# Reference dataset imputation helpers.
#
# Goal:
#   Use the empirical joint distribution in model-ready TMC controls to impute
#   coded missing values in population reference datasets. NFHS is used for the
#   premenopausal iCARE reference distribution and LASI is used for the
#   postmenopausal iCARE reference distribution.
#
# Current method:
#   External-donor random nearest-neighbor hot deck using StatMatch's Gower
#   distance machinery. TMC controls are the only donor pool; NFHS/LASI rows are
#   recipients. For each missingness pattern, observed covariates are used to
#   identify the 5 nearest TMC-control donors once per unique observed recipient
#   profile, then one donor is sampled at random and all missing covariates for
#   that recipient are copied jointly.
#
# Survey weights:
#   Each reference row also carries a `weight` column (NFHS DHS weight v005/1e6;
#   LASI all-India individual weight). The hot-deck step does not use it; it is
#   a passive passthrough so the iCARE step can supply it as
#   `model.ref.dataset.weights`, making the reference distribution population
#   representative.

standardize_missing <- function(x, missing_values = c(".", "", "NA", "N/A")) {
  if (inherits(x, "haven_labelled")) {
    x <- haven::zap_labels(x)
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    x[trimws(x) %in% missing_values] <- NA
  }

  x
}

numeric_code <- function(x) {
  suppressWarnings(as.numeric(as.character(standardize_missing(x))))
}

collapse_at_four <- function(x) {
  x <- numeric_code(x)
  ifelse(!is.na(x) & x >= 4, 4, x)
}

coerce_to_allowed_levels <- function(data, allowed_levels) {
  out <- data

  for (var in intersect(names(allowed_levels), names(out))) {
    out[[var]] <- numeric_code(out[[var]])
    out[[var]][!(out[[var]] %in% allowed_levels[[var]])] <- NA_real_
  }

  out
}

format_allowed_values <- function(constraints, covariates) {
  vapply(covariates, function(var) {
    if (!var %in% names(constraints)) {
      return(NA_character_)
    }

    paste(constraints[[var]], collapse = "/")
  }, character(1))
}

make_missing_report <- function(before_data, after_data, covariates, dataset,
                                constraints = list()) {
  before <- colSums(is.na(before_data[, covariates, drop = FALSE]))
  after <- colSums(is.na(after_data[, covariates, drop = FALSE]))

  data.frame(
    dataset = dataset,
    variable = covariates,
    allowed_imputed_values = format_allowed_values(constraints, covariates),
    missing_before = as.integer(before[covariates]),
    missing_after = as.integer(after[covariates]),
    missing_pct_before = round(100 * before[covariates] / nrow(before_data), 2),
    missing_pct_after = round(100 * after[covariates] / nrow(after_data), 2),
    row.names = NULL
  )
}

prepare_hotdeck_match_data <- function(data, vars) {
  allowed_levels <- reference_allowed_levels()
  categorical_vars <- intersect(names(allowed_levels), vars)
  match_data <- data[, vars, drop = FALSE]

  for (var in vars) {
    if (var %in% categorical_vars) {
      match_data[[var]] <- factor(numeric_code(match_data[[var]]),
                                  levels = allowed_levels[[var]])
    } else {
      match_data[[var]] <- numeric_code(match_data[[var]])
    }
  }

  match_data
}

compute_hotdeck_distances <- function(rec_match, don_match, dist_fun = "Gower") {
  if (tolower(dist_fun) != "gower") {
    stop("Only Gower distance is currently supported for cached hot-deck candidates.")
  }

  distances <- StatMatch::gower.dist(data.x = rec_match, data.y = don_match)
  distances[is.nan(distances)] <- 1
  distances[is.na(distances)] <- 1
  distances
}

nearest_candidate_matrix <- function(distance_matrix, donor_ids, nearest_donors) {
  n_recipients <- nrow(distance_matrix)
  k <- min(nearest_donors, length(donor_ids))
  candidate_ids <- matrix(NA_integer_, nrow = n_recipients, ncol = k)
  candidate_distances <- matrix(NA_real_, nrow = n_recipients, ncol = k)

  work <- as.matrix(distance_matrix)
  work[!is.finite(work)] <- Inf
  rows_with_finite_donor <- rowSums(is.finite(work)) > 0

  if (any(!rows_with_finite_donor)) {
    stop("No finite donor distances were available for at least one recipient row.")
  }

  row_index <- seq_len(n_recipients)
  for (candidate_i in seq_len(k)) {
    available <- rowSums(is.finite(work)) > 0
    available_rows <- row_index[available]
    min_columns <- max.col(-work[available, , drop = FALSE], ties.method = "first")

    candidate_ids[available_rows, candidate_i] <- donor_ids[min_columns]
    candidate_distances[available_rows, candidate_i] <-
      work[cbind(available_rows, min_columns)]
    work[cbind(available_rows, min_columns)] <- Inf
  }

  list(ids = candidate_ids, distances = candidate_distances)
}

build_hotdeck_candidates <- function(reference_data,
                                     donor_data,
                                     covariates,
                                     match_vars,
                                     constrained_values = list(),
                                     seed = 20260507,
                                     nearest_donors = 5,
                                     chunk_size = 5000,
                                     dist_fun = "Gower",
                                     verbose = FALSE) {
  if (!requireNamespace("StatMatch", quietly = TRUE)) {
    stop(
      "The StatMatch package is required for random nearest-neighbor hot-deck imputation. ",
      "Install it with install.packages('StatMatch')."
    )
  }

  stopifnot(all(covariates %in% names(reference_data)))
  stopifnot(all(covariates %in% names(donor_data)))

  set.seed(seed)
  out <- reference_data
  match_vars <- unique(match_vars)
  match_vars <- intersect(match_vars, intersect(names(reference_data), names(donor_data)))
  donor_required_vars <- unique(c(covariates, match_vars))
  donor_complete <- donor_data[
    complete.cases(donor_data[, donor_required_vars, drop = FALSE]),
    donor_required_vars,
    drop = FALSE
  ]

  if (nrow(donor_complete) == 0) {
    stop("No complete donor rows are available for hot-deck imputation.")
  }

  donor_complete$.donor_id <- seq_len(nrow(donor_complete))
  candidate_groups <- list()
  candidate_diagnostics <- list()

  make_candidate_group <- function(recipient_rows, donor_pool, active_match_vars,
                                   signature_label) {
    if (length(active_match_vars) > 0) {
      donor_pool <- donor_pool[
        complete.cases(donor_pool[, active_match_vars, drop = FALSE]),
        ,
        drop = FALSE
      ]
    }

    if (nrow(donor_pool) == 0) {
      stop("No donors available for missing signature: ", signature_label)
    }

    k_nearest <- min(nearest_donors, nrow(donor_pool))
    candidate_ids <- matrix(NA_integer_, nrow = length(recipient_rows), ncol = k_nearest)
    candidate_distances <- matrix(NA_real_, nrow = length(recipient_rows), ncol = k_nearest)

    if (length(active_match_vars) == 0) {
      chunks <- split(seq_along(recipient_rows),
                      ceiling(seq_along(recipient_rows) / chunk_size))

      for (chunk_i in seq_along(chunks)) {
        chunk_positions <- chunks[[chunk_i]]
        set.seed(seed + chunk_i)
        sampled_pool_rows <- matrix(
          unlist(replicate(
            length(chunk_positions),
            sample(seq_len(nrow(donor_pool)), k_nearest, replace = FALSE),
            simplify = FALSE
          )),
          nrow = length(chunk_positions),
          ncol = k_nearest,
          byrow = TRUE
        )
        candidate_ids[chunk_positions, ] <- matrix(
          donor_pool$.donor_id[as.vector(sampled_pool_rows)],
          nrow = length(chunk_positions),
          ncol = k_nearest
        )
        candidate_distances[chunk_positions, ] <- NA_real_
      }

      n_unique_profiles <- NA_integer_
    } else {
      rec_match_all <- prepare_hotdeck_match_data(
        reference_data[recipient_rows, , drop = FALSE],
        active_match_vars
      )
      profile_frame <- data.frame(
        lapply(rec_match_all, as.character),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      profile_keys <- do.call(paste, c(profile_frame, sep = "\r"))
      unique_profile_rows <- !duplicated(profile_keys)
      unique_profile_keys <- profile_keys[unique_profile_rows]
      profile_index <- match(profile_keys, unique_profile_keys)
      unique_match <- rec_match_all[unique_profile_rows, , drop = FALSE]

      unique_candidate_ids <- matrix(NA_integer_,
                                     nrow = nrow(unique_match),
                                     ncol = k_nearest)
      unique_candidate_distances <- matrix(NA_real_,
                                           nrow = nrow(unique_match),
                                           ncol = k_nearest)
      profile_chunks <- split(seq_len(nrow(unique_match)),
                              ceiling(seq_len(nrow(unique_match)) / chunk_size))
      don_match <- prepare_hotdeck_match_data(donor_pool, active_match_vars)

      for (chunk_i in seq_along(profile_chunks)) {
        chunk_profiles <- profile_chunks[[chunk_i]]
        distances <- compute_hotdeck_distances(
          unique_match[chunk_profiles, , drop = FALSE],
          don_match,
          dist_fun = dist_fun
        )
        candidates <- nearest_candidate_matrix(
          distances,
          donor_ids = donor_pool$.donor_id,
          nearest_donors = k_nearest
        )
        unique_candidate_ids[chunk_profiles, ] <- candidates$ids
        unique_candidate_distances[chunk_profiles, ] <- candidates$distances
      }

      candidate_ids <- unique_candidate_ids[profile_index, , drop = FALSE]
      candidate_distances <- unique_candidate_distances[profile_index, , drop = FALSE]
      n_unique_profiles <- nrow(unique_match)
    }

    summarize_distance <- function(x, fun) {
      if (all(is.na(x))) {
        return(NA_real_)
      }
      fun(x, na.rm = TRUE)
    }

    candidate_diagnostics[[length(candidate_diagnostics) + 1]] <<- data.frame(
      missing_signature = signature_label,
      active_match_vars = paste(active_match_vars, collapse = ";"),
      n_recipients = length(recipient_rows),
      n_unique_profiles = n_unique_profiles,
      n_donors_available = nrow(donor_pool),
      nearest_donors = k_nearest,
      mean_nearest_distance = summarize_distance(candidate_distances[, 1], mean),
      median_nearest_distance = summarize_distance(candidate_distances[, 1], median),
      mean_kth_distance = summarize_distance(candidate_distances[, k_nearest], mean),
      row.names = NULL
    )

    list(
      missing_signature = signature_label,
      missing_vars = strsplit(signature_label, "\t", fixed = TRUE)[[1]],
      recipient_rows = recipient_rows,
      candidate_ids = candidate_ids,
      candidate_distances = candidate_distances,
      active_match_vars = active_match_vars,
      n_donors_available = nrow(donor_pool),
      nearest_donors = k_nearest
    )
  }

  remaining <- which(!complete.cases(out[, covariates, drop = FALSE]))
  if (length(remaining) > 0) {
    missing_signatures <- vapply(remaining, function(row_i) {
      paste(covariates[is.na(out[row_i, covariates])], collapse = "\t")
    }, character(1))

    for (signature_i in seq_along(unique(missing_signatures))) {
      signature <- unique(missing_signatures)[signature_i]
      rows_for_signature <- remaining[missing_signatures == signature]
      missing_vars <- strsplit(signature, "\t", fixed = TRUE)[[1]]
      missing_vars <- missing_vars[nzchar(missing_vars)]

      if (verbose) {
        message(
          "Building hot-deck candidates for missing signature: ",
          gsub("\t", " + ", signature),
          " (n = ", length(rows_for_signature), ")"
        )
      }

      donor_candidates <- donor_complete
      for (var in intersect(missing_vars, names(constrained_values))) {
        donor_candidates <- donor_candidates[
          donor_candidates[[var]] %in% constrained_values[[var]],
          ,
          drop = FALSE
        ]
      }

      active_match_vars <- setdiff(match_vars, missing_vars)
      active_match_vars <- active_match_vars[vapply(active_match_vars, function(var) {
        all(!is.na(out[[var]][rows_for_signature]))
      }, logical(1))]

      candidate_groups[[length(candidate_groups) + 1]] <- make_candidate_group(
        recipient_rows = rows_for_signature,
        donor_pool = donor_candidates,
        active_match_vars = active_match_vars,
        signature_label = signature
      )
    }
  }

  diagnostics <- if (length(candidate_diagnostics) > 0) {
    do.call(rbind, candidate_diagnostics)
  } else {
    data.frame(
      missing_signature = character(0),
      active_match_vars = character(0),
      n_recipients = integer(0),
      n_unique_profiles = integer(0),
      n_donors_available = integer(0),
      nearest_donors = integer(0),
      mean_nearest_distance = numeric(0),
      median_nearest_distance = numeric(0),
      mean_kth_distance = numeric(0)
    )
  }

  list(
    reference_data = reference_data,
    donor_complete = donor_complete,
    covariates = covariates,
    match_vars = match_vars,
    constrained_values = constrained_values,
    nearest_donors = nearest_donors,
    dist_fun = dist_fun,
    candidate_groups = candidate_groups,
    diagnostics = diagnostics
  )
}

impute_from_hotdeck_candidates <- function(candidate_object, seed = 20260507) {
  set.seed(seed)
  out <- candidate_object$reference_data
  donor_complete <- candidate_object$donor_complete

  for (group in candidate_object$candidate_groups) {
    candidate_ids <- group$candidate_ids
    candidate_counts <- rowSums(!is.na(candidate_ids))

    if (any(candidate_counts == 0)) {
      stop("At least one recipient has no hot-deck donor candidates.")
    }

    sampled_columns <- vapply(candidate_counts, function(n_candidates) {
      sample.int(n_candidates, 1)
    }, integer(1))
    sampled_donor_ids <- candidate_ids[cbind(seq_len(nrow(candidate_ids)), sampled_columns)]

    for (var in group$missing_vars) {
      out[[var]][group$recipient_rows] <- donor_complete[[var]][sampled_donor_ids]
    }
  }

  remaining <- which(!complete.cases(out[, candidate_object$covariates, drop = FALSE]))
  if (length(remaining) > 0) {
    stop("Hot-deck imputation left missing rows after donor-candidate sampling.")
  }

  attr(out, "hotdeck_diagnostics") <- candidate_object$diagnostics
  out
}

impute_joint_hotdeck <- function(reference_data,
                                 donor_data,
                                 covariates,
                                 match_vars,
                                 constrained_values = list(),
                                 seed = 20260507,
                                 nearest_donors = 5,
                                 chunk_size = 5000,
                                 dist_fun = "Gower",
                                 verbose = FALSE) {
  candidates <- build_hotdeck_candidates(
    reference_data = reference_data,
    donor_data = donor_data,
    covariates = covariates,
    match_vars = match_vars,
    constrained_values = constrained_values,
    seed = seed,
    nearest_donors = nearest_donors,
    chunk_size = chunk_size,
    dist_fun = dist_fun,
    verbose = verbose
  )

  impute_from_hotdeck_candidates(candidates, seed = seed)
}

prepare_control_donors <- function(control_data) {
  donors <- data.frame(
    menopause = as.character(control_data$menopause),
    Age = numeric_code(control_data$Age),
    bmicat = numeric_code(control_data$bmicat),
    catht = numeric_code(control_data$catht),
    ratio = numeric_code(control_data$ratio),
    ageatfirstfulltermpreg_cat = numeric_code(control_data$ageatfirstfulltermpreg_cat),
    fulltermpreg_cat = numeric_code(control_data$fulltermpreg_cat),
    totalmiscarriage_cat = numeric_code(control_data$totalmiscarriage_cat),
    residancefirsttwentyyears = numeric_code(control_data$residancefirsttwentyyears),
    br_lump_yn = numeric_code(control_data$br_lump_yn),
    tobacco_chewing_yn = numeric_code(control_data$tobacco_chewing_yn),
    stringsAsFactors = FALSE
  )

  coerce_to_allowed_levels(donors, reference_allowed_levels())
}

prepare_nfhs_reference <- function(path) {
  nfhs <- haven::read_dta(path)
  nfhs[] <- lapply(nfhs, standardize_missing)

  # Column names follow the 2026-05 Sneha NFHS extract (NFHS_ref_main.dta).
  # Covariates are already coded to the model categories. `v005` is the DHS
  # women's individual sample weight; it is rescaled by 1e6 per DHS
  # convention so the stored weight is human-readable (the iCARE weighted
  # mean is scale-invariant, so the rescaling does not change any result).
  out <- data.frame(
    id = nfhs$caseid,
    weight = numeric_code(nfhs$v005) / 1e6,
    Age = numeric_code(nfhs$v012),
    bmicat = numeric_code(nfhs$bmi_cat),
    catht = numeric_code(nfhs$height_cat),
    ratio = numeric_code(nfhs$whr_cat),
    ageatfirstfulltermpreg_cat = numeric_code(nfhs$age_at_fftp),
    totalmiscarriage_cat = numeric_code(nfhs$miscarriage),
    residancefirsttwentyyears = numeric_code(nfhs$residence),
    br_lump_yn = numeric_code(nfhs$breast_lump),
    tobacco_chewing_yn = numeric_code(nfhs$tobacco_chewing),
    stringsAsFactors = FALSE
  )

  coerce_to_allowed_levels(out, reference_allowed_levels())
}

prepare_lasi_reference <- function(path, age_path) {
  lasi <- haven::read_dta(path)
  lasi[] <- lapply(lasi, standardize_missing)

  # Column names follow the 2026-05 Sneha LASI extract (lasi_new.dta).
  # Covariates are already coded to the model categories and
  # `indiaindividualweight` is the LASI all-India individual weight.
  #
  # The Sneha LASI extract does not carry respondent age, so age is merged
  # from the prior LASI extract (LASI_extracted2.dta, variable `dm005`) by
  # `prim_key`. The two files share `prim_key` exactly.
  age_source <- haven::read_dta(age_path)
  age_source[] <- lapply(age_source, standardize_missing)
  age_lookup <- data.frame(
    prim_key = as.character(age_source$prim_key),
    Age = numeric_code(age_source$dm005),
    stringsAsFactors = FALSE
  )
  age_lookup <- age_lookup[!duplicated(age_lookup$prim_key), , drop = FALSE]
  matched_age <- age_lookup$Age[match(as.character(lasi$prim_key),
                                      age_lookup$prim_key)]

  if (anyNA(matched_age)) {
    warning(sum(is.na(matched_age)),
            " LASI rows had no prim_key match for age in ", age_path)
  }

  out <- data.frame(
    id = lasi$prim_key,
    weight = numeric_code(lasi$indiaindividualweight),
    Age = matched_age,
    bmicat = numeric_code(lasi$bmi_cat),
    catht = numeric_code(lasi$height_cat),
    ratio = numeric_code(lasi$whr_cat),
    fulltermpreg_cat = collapse_at_four(lasi$full_term_preg),
    totalmiscarriage_cat = numeric_code(lasi$miscarriages),
    residancefirsttwentyyears = numeric_code(lasi$residence_hist),
    br_lump_yn = numeric_code(lasi$breast_lump),
    tobacco_chewing_yn = numeric_code(lasi$tob_chew),
    stringsAsFactors = FALSE
  )

  coerce_to_allowed_levels(out, reference_allowed_levels())
}

reference_allowed_levels <- function() {
  list(
    bmicat = c(1, 2, 3),
    catht = c(0, 1, 2, 3),
    ratio = c(0, 1, 2),
    ageatfirstfulltermpreg_cat = c(1, 2, 3, 4, 5, 6),
    fulltermpreg_cat = c(0, 1, 2, 3, 4),
    totalmiscarriage_cat = c(0, 1, 2),
    residancefirsttwentyyears = c(1, 2),
    br_lump_yn = c(0, 1),
    tobacco_chewing_yn = c(0, 1, 2)
  )
}

reference_missing_constraints <- function(dataset) {
  allowed <- reference_allowed_levels()

  if (identical(dataset, "NFHS_premenopausal")) {
    allowed$totalmiscarriage_cat <- c(1, 2)
    allowed$tobacco_chewing_yn <- c(1, 2)
    return(allowed)
  }

  if (identical(dataset, "LASI_postmenopausal")) {
    allowed$tobacco_chewing_yn <- c(1, 2)
    return(allowed)
  }

  stop("Unsupported reference dataset for missingness constraints: ", dataset)
}

premeno_reference_covariates <- function() {
  c(
    "bmicat", "ratio", "residancefirsttwentyyears",
    "ageatfirstfulltermpreg_cat", "catht", "totalmiscarriage_cat",
    "br_lump_yn", "tobacco_chewing_yn"
  )
}

postmeno_reference_covariates <- function() {
  c(
    "bmicat", "ratio", "residancefirsttwentyyears",
    "fulltermpreg_cat", "catht", "totalmiscarriage_cat",
    "br_lump_yn", "tobacco_chewing_yn"
  )
}

premeno_match_vars <- function() {
  c("Age", premeno_reference_covariates())
}

postmeno_match_vars <- function() {
  c("Age", postmeno_reference_covariates())
}

make_level_frequencies <- function(data, covariates, allowed_levels) {
  do.call(rbind, lapply(covariates, function(var) {
    levels <- allowed_levels[[var]]
    x <- factor(numeric_code(data[[var]]), levels = levels)
    counts <- table(x, useNA = "no")
    total <- sum(counts)

    data.frame(
      variable = var,
      level = names(counts),
      n = as.integer(counts),
      pct = if (total > 0) round(100 * as.numeric(counts) / total, 2) else NA_real_,
      total = as.integer(total),
      row.names = NULL
    )
  }))
}

make_reference_frequency_diagnostics <- function(reference_outputs,
                                                 control_data,
                                                 output_dir = NULL) {
  allowed <- reference_allowed_levels()
  control_donors <- prepare_control_donors(control_data)

  pre_covariates <- premeno_reference_covariates()
  post_covariates <- postmeno_reference_covariates()

  nfhs_premeno <- reference_outputs$nfhs_premeno[
    reference_outputs$nfhs_premeno$Age >= 25 &
      reference_outputs$nfhs_premeno$Age < 45,
    ,
    drop = FALSE
  ]
  lasi_postmeno <- reference_outputs$lasi_postmeno[
    reference_outputs$lasi_postmeno$Age >= 45 &
      reference_outputs$lasi_postmeno$Age <= 68,
    ,
    drop = FALSE
  ]

  pre_donors <- control_donors[
    control_donors$menopause == "premenopausal" &
      complete.cases(control_donors[, pre_covariates, drop = FALSE]),
    ,
    drop = FALSE
  ]
  post_donors <- control_donors[
    control_donors$menopause == "postmenopausal" &
      complete.cases(control_donors[, post_covariates, drop = FALSE]),
    ,
    drop = FALSE
  ]

  compare_group <- function(reference_data, donor_data, covariates, model_group,
                            reference_dataset, donor_dataset) {
    reference_freq <- make_level_frequencies(reference_data, covariates, allowed)
    donor_freq <- make_level_frequencies(donor_data, covariates, allowed)

    names(reference_freq)[names(reference_freq) %in% c("n", "pct", "total")] <-
      c("reference_n", "reference_pct", "reference_total")
    names(donor_freq)[names(donor_freq) %in% c("n", "pct", "total")] <-
      c("donor_n", "donor_pct", "donor_total")

    out <- merge(reference_freq, donor_freq, by = c("variable", "level"),
                 all = TRUE, sort = FALSE)
    out$model_group <- model_group
    out$reference_dataset <- reference_dataset
    out$donor_dataset <- donor_dataset
    out$pct_difference <- round(out$reference_pct - out$donor_pct, 2)
    out$abs_pct_difference <- abs(out$pct_difference)

    out[, c(
      "model_group", "reference_dataset", "donor_dataset",
      "variable", "level",
      "reference_n", "reference_pct", "reference_total",
      "donor_n", "donor_pct", "donor_total",
      "pct_difference", "abs_pct_difference"
    )]
  }

  diagnostics <- rbind(
    compare_group(
      nfhs_premeno, pre_donors, pre_covariates,
      model_group = "premenopausal",
      reference_dataset = "NFHS_premenopausal_imputed_age_25_44",
      donor_dataset = "TMC_premenopausal_complete_control_donors"
    ),
    compare_group(
      lasi_postmeno, post_donors, post_covariates,
      model_group = "postmenopausal",
      reference_dataset = "LASI_postmenopausal_imputed_age_45_68",
      donor_dataset = "TMC_postmenopausal_complete_control_donors"
    )
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(
      diagnostics,
      file.path(output_dir, "reference_imputation_frequency_diagnostics.csv"),
      row.names = FALSE
    )
  }

  diagnostics
}

build_reference_datasets <- function(control_data,
                                     nfhs_path,
                                     lasi_path,
                                     lasi_age_path,
                                     output_dir = NULL,
                                     seed = 20260507,
                                     nearest_donors = 5,
                                     chunk_size = 5000,
                                     m = 1,
                                     verbose = FALSE) {
  if (m < 1) {
    stop("m must be at least 1.")
  }

  control_donors <- prepare_control_donors(control_data)
  pre_donors <- control_donors[control_donors$menopause == "premenopausal", ]
  post_donors <- control_donors[control_donors$menopause == "postmenopausal", ]

  nfhs <- prepare_nfhs_reference(nfhs_path)
  lasi <- prepare_lasi_reference(lasi_path, lasi_age_path)

  pre_covariates <- premeno_reference_covariates()
  post_covariates <- postmeno_reference_covariates()

  nfhs_before <- nfhs
  lasi_before <- lasi
  nfhs_constraints <- reference_missing_constraints("NFHS_premenopausal")
  lasi_constraints <- reference_missing_constraints("LASI_postmenopausal")

  nfhs_candidates <- build_hotdeck_candidates(
    reference_data = nfhs,
    donor_data = pre_donors,
    covariates = pre_covariates,
    match_vars = premeno_match_vars(),
    constrained_values = nfhs_constraints,
    seed = seed,
    nearest_donors = nearest_donors,
    chunk_size = chunk_size,
    verbose = verbose
  )

  lasi_candidates <- build_hotdeck_candidates(
    reference_data = lasi,
    donor_data = post_donors,
    covariates = post_covariates,
    match_vars = postmeno_match_vars(),
    constrained_values = lasi_constraints,
    seed = seed + 1,
    nearest_donors = nearest_donors,
    chunk_size = chunk_size,
    verbose = verbose
  )

  build_panel <- function(panel_i) {
    panel_seed <- seed + (panel_i - 1) * 1000
    nfhs_panel <- impute_from_hotdeck_candidates(nfhs_candidates, seed = panel_seed)
    lasi_panel <- impute_from_hotdeck_candidates(lasi_candidates, seed = panel_seed + 1)

    diagnostics <- rbind(
      cbind(
        panel = panel_i,
        make_missing_report(
          nfhs_before, nfhs_panel, pre_covariates, "NFHS_premenopausal",
          constraints = nfhs_constraints
        )
      ),
      cbind(
        panel = panel_i,
        make_missing_report(
          lasi_before, lasi_panel, post_covariates, "LASI_postmenopausal",
          constraints = lasi_constraints
        )
      )
    )

    list(
      panel = panel_i,
      nfhs_premeno = nfhs_panel[, c("id", "Age", "weight", pre_covariates),
                                drop = FALSE],
      lasi_postmeno = lasi_panel[, c("id", "Age", "weight", post_covariates),
                                 drop = FALSE],
      diagnostics = diagnostics
    )
  }

  reference_panels <- lapply(seq_len(m), build_panel)
  diagnostics <- do.call(rbind, lapply(reference_panels, `[[`, "diagnostics"))

  donor_counts <- data.frame(
    donor_group = c("premenopausal_controls", "postmenopausal_controls"),
    n = c(nrow(pre_donors), nrow(post_donors)),
    complete_n = c(
      sum(complete.cases(pre_donors[, pre_covariates, drop = FALSE])),
      sum(complete.cases(post_donors[, post_covariates, drop = FALSE]))
    ),
    method = "Cached random nearest-neighbor hot deck using StatMatch::gower.dist",
    nearest_donors = nearest_donors,
    imputations = m
  )

  hotdeck_diagnostics <- rbind(
    cbind(dataset = "NFHS_premenopausal", nfhs_candidates$diagnostics),
    cbind(dataset = "LASI_postmenopausal", lasi_candidates$diagnostics)
  )

  frequency_diagnostics <- do.call(rbind, lapply(reference_panels, function(panel) {
    panel_result <- list(
      nfhs_premeno = panel$nfhs_premeno,
      lasi_postmeno = panel$lasi_postmeno
    )
    cbind(
      panel = panel$panel,
      make_reference_frequency_diagnostics(
        panel_result,
        control_data = control_data
      )
    )
  }))

  result <- list(
    nfhs_premeno = reference_panels[[1]]$nfhs_premeno,
    lasi_postmeno = reference_panels[[1]]$lasi_postmeno,
    reference_panels = reference_panels,
    diagnostics = diagnostics,
    donor_counts = donor_counts,
    hotdeck_diagnostics = hotdeck_diagnostics,
    frequency_diagnostics = frequency_diagnostics,
    candidate_metadata = list(
      nearest_donors = nearest_donors,
      chunk_size = chunk_size,
      seed = seed,
      imputations = m,
      nfhs_candidate_groups = length(nfhs_candidates$candidate_groups),
      lasi_candidate_groups = length(lasi_candidates$candidate_groups)
    )
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(result$nfhs_premeno, file.path(output_dir, "nfhs_premeno_reference_imputed.rds"))
    saveRDS(result$lasi_postmeno, file.path(output_dir, "lasi_postmeno_reference_imputed.rds"))
    saveRDS(reference_panels, file.path(output_dir, "reference_imputation_panels.rds"))
    saveRDS(
      list(nfhs = nfhs_candidates, lasi = lasi_candidates),
      file.path(output_dir, "reference_hotdeck_candidates.rds")
    )

    for (panel in reference_panels) {
      panel_suffix <- sprintf("panel_%02d", panel$panel)
      saveRDS(
        panel$nfhs_premeno,
        file.path(output_dir, paste0("nfhs_premeno_reference_imputed_", panel_suffix, ".rds"))
      )
      saveRDS(
        panel$lasi_postmeno,
        file.path(output_dir, paste0("lasi_postmeno_reference_imputed_", panel_suffix, ".rds"))
      )
    }

    saveRDS(result, file.path(output_dir, "reference_imputation_outputs.rds"))
    write.csv(diagnostics, file.path(output_dir, "reference_imputation_missingness.csv"),
              row.names = FALSE)
    write.csv(donor_counts, file.path(output_dir, "reference_imputation_donors.csv"),
              row.names = FALSE)
    write.csv(hotdeck_diagnostics,
              file.path(output_dir, "reference_imputation_hotdeck_diagnostics.csv"),
              row.names = FALSE)
    write.csv(frequency_diagnostics,
              file.path(output_dir, "reference_imputation_frequency_diagnostics.csv"),
              row.names = FALSE)
  }

  result
}

# Probabilistic hot-deck imputation, stratified by `group_vars`.
#
# Within each stratum, donors are the rows complete in `covariates` and
# recipients are rows with any missing value. The hot deck is run `m` times
# (5 by default) using the same Gower-distance random nearest-neighbor engine
# as `impute_joint_hotdeck`. Imputed values for each originally-missing cell
# are pooled across the `m` panels: modal category for categorical covariates,
# row-mean for any covariates listed in `continuous_vars`.
#
# Used by breast_cancer_risk_model.Rmd for both the TMC training-data imputation
# (stratified by case status x menopause) and the validation-data imputation
# (stratified by menopause only, since case/control status is not available
# at prediction time).
impute_joint_hotdeck_grouped <- function(data,
                                          covariates,
                                          group_vars,
                                          continuous_vars = c("Age"),
                                          nearest_donors = 5,
                                          m = 5,
                                          seed = 20260529,
                                          verbose = FALSE) {
  stopifnot(all(group_vars %in% names(data)))
  stopifnot(all(covariates %in% names(data)))

  data <- as.data.frame(data)
  group_keys <- do.call(paste, c(data[, group_vars, drop = FALSE], sep = "|"))
  unique_groups <- unique(group_keys[!is.na(group_keys)])

  imputed_panels <- vector("list", m)

  for (panel_i in seq_len(m)) {
    panel_data <- data

    for (g_idx in seq_along(unique_groups)) {
      g <- unique_groups[g_idx]
      rows <- which(group_keys == g)
      if (length(rows) == 0) next
      sub <- data[rows, , drop = FALSE]

      donor_complete <- sub[complete.cases(sub[, covariates, drop = FALSE]),
                            , drop = FALSE]
      has_missing <- !complete.cases(sub[, covariates, drop = FALSE])

      if (!any(has_missing)) next
      if (nrow(donor_complete) == 0) {
        if (verbose) {
          message("No complete donors in stratum '", g, "'; skipping ",
                  sum(has_missing), " recipients")
        }
        next
      }

      sub_imputed <- impute_joint_hotdeck(
        reference_data = sub,
        donor_data = donor_complete,
        covariates = covariates,
        match_vars = unique(c(continuous_vars, covariates)),
        seed = seed + panel_i * 1000L + g_idx,
        nearest_donors = nearest_donors,
        verbose = verbose
      )
      panel_data[rows, ] <- sub_imputed
    }

    imputed_panels[[panel_i]] <- panel_data
  }

  pooled <- data
  for (var in covariates) {
    was_missing <- is.na(data[[var]])
    if (!any(was_missing)) next

    if (var %in% continuous_vars && is.numeric(data[[var]])) {
      panel_vals <- vapply(imputed_panels,
                            function(p) p[[var]][was_missing],
                            numeric(sum(was_missing)))
      if (sum(was_missing) == 1L) {
        pooled[[var]][was_missing] <- mean(panel_vals, na.rm = TRUE)
      } else {
        pooled[[var]][was_missing] <- rowMeans(panel_vals, na.rm = TRUE)
      }
    } else {
      panel_chars <- vapply(imputed_panels,
                             function(p) as.character(p[[var]][was_missing]),
                             character(sum(was_missing)))
      panel_chars <- matrix(panel_chars,
                             nrow = sum(was_missing),
                             ncol = m)
      modal_chr <- apply(panel_chars, 1, function(row) {
        non_na <- row[!is.na(row)]
        if (length(non_na) == 0) return(NA_character_)
        names(sort(table(non_na), decreasing = TRUE))[1]
      })
      new_vals <- if (is.numeric(data[[var]])) {
        suppressWarnings(as.numeric(modal_chr))
      } else if (is.factor(data[[var]])) {
        factor(modal_chr, levels = levels(data[[var]]))
      } else {
        modal_chr
      }
      pooled[[var]][was_missing] <- new_vals
    }
  }

  # Fallback pass: a stratum may have had zero complete-case donors (e.g. when
  # the group key itself is missing), leaving a small number of rows still
  # incomplete. Mirror the legacy median/mode helper's overall-distribution
  # fallback by running one more hot-deck against the global complete-case
  # pool for any rows still carrying NAs.
  still_missing <- which(!complete.cases(pooled[, covariates, drop = FALSE]))
  if (length(still_missing) > 0) {
    overall_donors <- pooled[
      complete.cases(pooled[, covariates, drop = FALSE]),
      , drop = FALSE
    ]
    if (nrow(overall_donors) > 0) {
      stragglers <- impute_joint_hotdeck(
        reference_data = pooled[still_missing, , drop = FALSE],
        donor_data = overall_donors,
        covariates = covariates,
        match_vars = unique(c(continuous_vars, covariates)),
        seed = seed + 9999L,
        nearest_donors = nearest_donors,
        verbose = verbose
      )
      pooled[still_missing, ] <- stragglers
      if (verbose) {
        message("Overall-distribution fallback imputed ",
                length(still_missing), " straggler rows")
      }
    } else if (verbose) {
      message(length(still_missing),
              " rows could not be imputed -- no overall donors available")
    }
  }

  attr(pooled, "hotdeck_panels") <- m
  attr(pooled, "hotdeck_group_vars") <- group_vars
  attr(pooled, "hotdeck_method") <-
    "random nearest-neighbor hot deck (Gower distance) with modal pooling + overall fallback"
  pooled
}
