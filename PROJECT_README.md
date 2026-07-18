# TMC Breast Cancer Absolute Risk Model — Project README

End‑to‑end documentation of the TMC breast cancer risk project: from the raw
datasets, through cohort assembly, covariate construction, model estimation and
reference building, to the calibrated **absolute risk** model and its outputs.

> Scope: the analysis pipeline in `TMC-BrCaProj/code/analysis/breast_cancer_risk_model.Rmd`
> (a.k.a. `Saket_code.Rmd`). The polygenic score (PRS) work is documented
> separately and is **excluded** here.

---

## Index

1. [Project overview](#1-project-overview)
2. [The datasets at hand](#2-the-datasets-at-hand)
3. [Cohort assembly — from files to the model frame](#3-cohort-assembly--from-files-to-the-model-frame)
4. [Covariate definitions and categorisation](#4-covariate-definitions-and-categorisation)
5. [Menopause stratification (age 48)](#5-menopause-stratification-age-48)
6. [Missing‑data imputation](#6-missing-data-imputation)
7. [Relative‑risk model estimation](#7-relative-risk-model-estimation)
8. [Reference population construction (NFHS / LASI)](#8-reference-population-construction-nfhs--lasi)
9. [Incidence and competing mortality](#9-incidence-and-competing-mortality)
10. [Absolute‑risk computation (iCARE)](#10-absolute-risk-computation-icare)
11. [Population projection](#11-population-projection)
12. [Validation](#12-validation)
13. [Outputs and the model object](#13-outputs-and-the-model-object)
14. [Cohort sizes (the funnel)](#14-cohort-sizes-the-funnel)
15. [Repository structure](#15-repository-structure)
16. [Known limitations and caveats](#16-known-limitations-and-caveats)
17. [How to run](#17-how-to-run)

---

## 1. Project overview

The goal is an **absolute (cumulative) breast cancer risk** model for Indian
women — the probability of developing breast cancer over 5 years, 10 years, and
lifetime (to age 80), calibrated to Indian incidence and reported with a
population percentile.

The statistical engine is **iCARE** (Individualized Coherent Absolute Risk
Estimator), which combines three ingredients:

1. **Relative risks** for a set of questionnaire covariates, estimated from the
   TMC case–control study (Section 7).
2. A **reference population** (NFHS‑5 for premenopausal, LASI for postmenopausal)
   that supplies the covariate distribution used to calibrate the model
   (Section 8).
3. **Age‑specific incidence and competing‑mortality rates** that anchor the
   relative risks to an absolute scale (Section 9).

The model is **stratified by menopausal status**: separate relative‑risk models
and separate references for premenopausal and postmenopausal women.

---

## 2. The datasets at hand

Every file below is read by the analysis. They fall into four functional groups.

### 2.1 Estimation cohort — TMC case–control data
The relative‑risk model is fit on the Tata Memorial Centre (TMH) case–control
study. Covariates arrive as separate per‑domain extracts merged on participant
`ID`, across two data vintages.

**Current drop — `BCRPM_Data_24.11` (Jan 2026), used in the merge:**

| File | Provides |
| --- | --- |
| `Anthro/Anthro_BCRPM_10.10.2025.xlsx` | height, weight, waist, hip, Age → BMI, WHR, height |
| `Personal_med/personal_med_BCRPM_11.10.2025.xlsx` | breast lump, diabetes, hypertension |
| `Chewing_incl/chewing_BCRPM_10.10.2025.xlsx` | tobacco chewing + age started (timing) |
| `Smoking_incl/Smoking_BCRPM_10.10.2025.xlsx` | smoking |
| `Alcohol/Alcohol_07.10.25.xlsx` | alcohol |

**Older extracts — `28.03.2025_Risk prediction model data` (2023–2024):**

| Status | Files |
| --- | --- |
| **Used** (no newer version) | `Reproductive/Reproductive_25.07.2023_final.xlsx` (parity, age at first pregnancy, OC, breastfeeding, abortions, miscarriages, menopause); `Demographics/Demography_details_…17062023.xlsx`; `Family_history_of_cancer/…` (2 files) |
| **Read but superseded** | 2023/2024 Alcohol (×3), Anthropometry (×2), Chewing, Personal_medical, Smoking; Residential_details set (`Last_5_3_year_RU`, `Last_5_3_year_state`, `POB_duration`, `Residential_birth_current`) |

**Baseline, residence, IDs, validation:**

| File | Role |
| --- | --- |
| `Questionnaire_28.04.2017.dta` | 2017 baseline cohort — rescues phase‑1 women (menopause/reproductive history) |
| `Article- …/Epi_data_residential/…/First_20_10_years_RU_9.11.25/Sheet1‑Table 1.csv` | first‑20‑years urban/rural residence |
| `Breast cancer data - CCE.TMC/data_new_meta_ids.xlsx` | "DOE" file: ID, Age, Phase, DateOfInterview, Study Centre, abortion/miscarriage counts → train/test split + TMH filter |
| `Breast cancer data - CCE.TMC/Meta_Data.xlsx` | Sample_ID crosswalk + GWAS/QC flags + Status |
| `Breast cancer data - CCE.TMC/mamogram data _24.5.2025.xlsx` | mammography‑linked IDs → date cutoff for the split |
| `Breast cancer data - CCE.TMC/Breast cancer _Validation data set.xlsx` | held‑out validation cohort |

### 2.2 Reference populations

| File | Role |
| --- | --- |
| `reference_data/Data-Sneha-reference/NFHS_ref_main.dta` | NFHS‑5 premenopausal reference |
| `reference_data/Data-Sneha-reference/lasi_new.dta` | LASI wave‑1 postmenopausal reference |
| `reference_data/LASI_extracted2.dta` | LASI age extract (merged for age) |

### 2.3 Rates and population

| File | Role |
| --- | --- |
| `BC_Incidence_Rates.xlsx` | age‑specific breast cancer incidence (31 Indian PBCRs), per 100,000 |
| `Overall_Mortality_Rates.xlsx` | competing (all‑cause) mortality, per 1,000 |
| `Indian population numbers by Age (Census-2011).xlsx` | census population, for projection |

### 2.4 External example profiles

| File | Role |
| --- | --- |
| `parichoy_br_data_6.2.23.xlsx` | Parichoy external dataset, used as illustrative iCARE test profiles |

---

## 3. Cohort assembly — from files to the model frame

The per‑domain extracts are harmonised and merged into a single analysis table:

1. **Column hygiene.** Each extract is cleaned (`janitor::clean_names`) and joined
   on `id`; duplicate columns across files are dropped so each variable appears once.
2. **Domain join.** The Jan‑2026 anthropometry, personal‑medical, chewing, smoking
   and alcohol extracts are `full_join`‑ed with the 2023 reproductive, demography
   and family‑history extracts to form `final_data`.
3. **Phase‑1 rescue via the 2017 questionnaire.** New‑cohort women whose 2025 files
   lack a matching 2023 reproductive row are recovered by `bind_rows` with the
   2017 questionnaire cohort, keeping the old record on overlap
   (`distinct(id, .keep_all = TRUE)`). This is what lifts coverage from the Jan‑2026
   drop alone to the full TMH cohort.
4. **Residence merge.** First‑20‑years urban/rural (`residancefirsttwentyyears`) is
   matched in from the residential‑history CSV (Nov 2025).
5. **New medical / lifestyle fields.** Breast lump (`br_lump_yn`) is matched from
   the personal‑medical extract; tobacco‑chewing timing is derived from the chewing
   extract (age started vs age at first pregnancy).
6. **Abortion/miscarriage counts** are pulled from the DOE meta file.
7. **Train/validation split.** Using `mamogram data`, the median interview date of
   mammography‑linked records defines a cutoff: women interviewed on/before it →
   training; after → validation. Analysis is restricted to **Study Centre == "TMH"**.
8. **Case‑exclusion for reverse causation.** Incident cases whose breast‑lump
   history dates within three years of interview are excluded from **both** training
   and validation (symmetric).

The result is the model frame `data`, carrying case status (`br`), menopause, and
all categorised covariates.

---

## 4. Covariate definitions and categorisation

A single set of helper functions defines every cutpoint (the "single source of
truth"), applied identically to the TMC data, the external profiles and the
NFHS/LASI reference.

| Covariate | Categories |
| --- | --- |
| BMI (`bmicat`) | 1 = <25; 2 = 25–29.9; 3 = ≥30 kg/m² |
| Waist‑hip ratio (`ratio`) | 0 = ≤0.84; 1 = 0.85–0.94; 2 = >0.94 |
| Height (`catht`) | 0 = ≤150; 1 = 151–155; 2 = 156–160; 3 = >160 cm |
| Age at first full‑term pregnancy (`ageatfirstfulltermpreg_cat`, **premeno**) | 1 = <20; 2 = 20–21; 3 = 22–23; 4 = 24–25; 5 = ≥26; 6 = nulliparous |
| Parity (`fulltermpreg_cat`, **postmeno**) | 0 = nulliparous; 1; **2 = reference**; 3; 4 = ≥4 |
| Miscarriages (`totalmiscarriage_cat`) | 0 = none; 1 = one; 2 = ≥2 |
| Residence, first 20 yrs (`residancefirsttwentyyears`) | **1 = urban; 2 = rural** |
| Tobacco chewing (`tobacco_chewing_yn`) | 0 = never; 1 = started before first pregnancy; 2 = started at/after |
| Breast lump (`br_lump_yn`) | 0 = no; 1 = yes |

`br_lump_yn` and `tobacco_chewing_yn` are the covariates added in the current
model (relative to the earlier version that used OC use, breastfeeding duration
and abortion count).

---

## 5. Menopause stratification (age 48)

Everything downstream splits into a premenopausal and a postmenopausal pipeline.

**Why stratify.** Breast‑cancer risk factors differ by menopausal status (e.g. BMI
is protective premenopausally, harmful postmenopausally), so the relative‑risk
models differ, and the reference populations differ (NFHS for premenopausal, LASI
for postmenopausal). A pooled model would be misspecified.

**Where menopausal status comes from.** It is *observed*, not inferred from age:
the TMC data carries a recorded menopause label (used for model fitting), and the
reference pools are already menopause‑pure at source — NFHS was filtered to
premenopausal women (menopausal excluded) and LASI to postmenopausal women. So the
age cutoff does **not** classify menopause in either the TMC data or the reference.

**Why an age cutoff (48) is nonetheless required — the age‑overlap / incidence
reason.** The two reference pools overlap in age: NFHS (premenopausal) runs up to
49 and LASI (postmenopausal) starts at 45, so ages **45–48 appear in both**. iCARE
(and the population projection) assigns the age‑specific incidence rate λ(t) to a
population at each age. If an age belongs to *both* strata, that single total λ(t)
would be applied to the premenopausal reference **and** the postmenopausal
reference at the same age — ambiguous for calibration and double‑counting in the
projection. Forcing the strata to be **age‑disjoint** removes the overlap:

- premenopausal reference / projection band: ages **25–48**
- postmenopausal reference / projection band: ages **49–68**

Now every age maps to exactly one population, so the incidence rate at each age is
assigned unambiguously (25–48 → premenopausal, 49–68 → postmenopausal). This is
also precisely what makes it valid to use a **single total incidence table** for
both strata (Section 9): no age is shared between them.

**Why 48 specifically.** It is the integer age that **minimises menopausal
misclassification** against the TMC women with a recorded status — 16.5% at 48 vs
19.5% at 45. The cutoff also allocates the age‑overlap band: raising it from 45 to
48 grew the NFHS premenopausal reference (377,200 → 431,790, gaining the
premenopausal 45–48‑year‑olds) and shrank the LASI postmenopausal reference
(27,816 → 21,834, dropping its 45–48‑year‑olds).

**What it is not.** The cutoff is not a covariate (there is no age term in the
relative‑risk model; age enters only through the incidence rates), and it is not
used for per‑individual scoring — a woman self‑reports menopausal status and is
scored against the already‑pure pre‑ or post‑menopausal reference.

---

## 6. Missing‑data imputation

Missing categorical covariates are completed by **probabilistic hot‑deck
imputation** (random nearest‑neighbour donors by Gower distance):

- **Training data:** `m = 5` imputed panels, stratified by case status × menopause,
  pooled to a completed dataset used to fit the GLM.
- **Reference data:** the same hot‑deck fills missing NFHS/LASI covariates, using
  TMC controls as external donors where required, constrained to allowed category
  levels (Section 8).

---

## 7. Relative‑risk model estimation

Two logistic regressions (one per stratum) on the completed TMC data. The fitted
GLM includes an `Age` term; the **iCARE** relative‑risk formula deliberately
excludes age.

**Premenopausal (iCARE formula):**

```
br ~ bmicat*ratio
   + residancefirsttwentyyears*bmicat
   + residancefirsttwentyyears*ratio
   + ageatfirstfulltermpreg_cat
   + catht + totalmiscarriage_cat + br_lump_yn + tobacco_chewing_yn
```

**Postmenopausal:** same, with `relevel(fulltermpreg_cat, ref = "2")` replacing
age at first pregnancy.

The fitted coefficients are mapped onto the iCARE design matrix
(`get_log_rr_for_icare`) to produce the log‑relative‑risk vectors
**`logRR_premeno`** (~26 terms) and **`logRR_postmeno`** (~25 terms), saved to
`Data/logRR_premeno.rds` / `Data/logRR_postmeno.rds`.

---

## 8. Reference population construction (NFHS / LASI)

The reference supplies the population covariate distribution used to calibrate
the baseline hazard and to place percentiles.

- **Sources:** NFHS‑5 (premenopausal), LASI wave‑1 (postmenopausal), built by
  `build_reference_datasets()` from the Sneha‑prepared `.dta` files (LASI age
  merged from `LASI_extracted2.dta`).
- **Harmonisation:** covariates categorised to match the model exactly; survey
  weights retained (NFHS `v005`); missingness filled by hot‑deck (5 panels), with
  TMC controls as donors.
- **Panels:** `make_icare_reference_panel()` selects `Age + weight + covariates`
  and filters to the stratum age band (25–48 premeno; 49–68 postmeno), producing
  **5 aligned imputed panels per stratum**, each a set of `{age, weight, covariates}`.
- **Sizes (after the age‑48 bands):** NFHS premeno **431,790** rows/panel; LASI
  postmeno **21,834** rows/panel.
- **Note:** breast lump and tobacco‑chewing timing are **not observed** in
  NFHS/LASI and are therefore imputed in the reference.

---

## 9. Incidence and competing mortality

- **Incidence:** `BC_Incidence_Rates.xlsx` — age‑specific breast cancer rates from
  31 Indian PBCRs, per 100,000. Formatted to per‑integer‑age and divided by 100,000.
- **Competing mortality:** `Overall_Mortality_Rates.xlsx` — per 1,000, formatted and
  divided by 1,000.
- A **single, total incidence table is used for both strata** (there is no
  menopause‑specific λ_pre/λ_post split; a mixture split was prototyped and
  disabled — see Section 16).

---

## 10. Absolute‑risk computation (iCARE)

For an individual with covariate profile X, iCARE computes absolute risk over
`[a, a+τ]` by:

1. Calibrating the baseline hazard so the model, averaged over the reference,
   reproduces population incidence: `λ₀(t) = λ(t) / E_ref[RR(X)]`.
2. Integrating the individual hazard `λ₀(t)·RR(X)` against survival from both
   breast cancer and competing mortality over the interval.

Because the reference is multiply imputed, absolute risk is computed **once per
reference panel and averaged across the 5 panels**
(`compute_absolute_risk_over_reference_panels`), for 5‑year, 10‑year and (in the
deployment) lifetime horizons.

---

## 11. Population projection

To estimate the number of future cases in the Indian population:

- Census‑2011 population is split into age bands (**premeno 25–48, postmeno 49–68**).
- Expected cases = population(band) × survey‑weighted **mean absolute risk** of the
  stratum.
- A smooth age‑at‑menopause mixture was considered but replaced by the hard 48/49
  band (Section 16).

Indicative outputs: ~1.52M expected 5‑year cases, ~3.36M 10‑year cases (grand
totals, stable across the cutoff change).

---

## 12. Validation

Internal validation on the held‑out TMC sample yields AUCs of **0.572
(premenopausal)** and **0.600 (postmenopausal)** — the low‑to‑moderate
discrimination expected of a questionnaire‑only risk model.

---

## 13. Outputs and the model object

- `Data/logRR_premeno.rds`, `Data/logRR_postmeno.rds` — relative‑risk vectors.
- Reference panels (`TMC-BrCaProj/data/reference/…`).
- Pooled mean absolute risks — premeno 5y ≈ 0.159%, 10y ≈ 0.397%; postmeno 5y ≈
  0.596%, 10y ≈ 1.204%.
- The **deployment bundle** `bc_risk_input.rds` (assembled for the app): a list of
  `ref_premeno` (5 NFHS panels), `ref_postmeno` (5 LASI panels), `bc_incidence`,
  `mortality`, `logrr_pre`, `logrr_post`.

---

## 14. Cohort sizes (the funnel)

| Step | Count |
| --- | --- |
| BCRPM enrolment (all centres) | 7,747 |
| TMH centre (geographic eligibility) | 6,201 |
| Final model frame (training + validation) | 6,125 |
| — Training | 4,564 (premeno 2,403; postmeno 2,161) |
| — Validation | 1,561 (premeno 735; postmeno 826) |

Only ~76 TMH women fall out (breast‑lump 3‑year exclusion + unrecorded menopause).

---

## 15. Repository structure

```
TMC_BCProject/
  Data/                         raw datasets, reference data, rate tables, model files
    reference_data/             NFHS / LASI reference .dta + variable guides
    BCRPM_Data_24.11/           Jan-2026 covariate extracts (current)
    Breast cancer data - CCE.TMC/  older extracts, meta/ID files, validation set
  Article- Breast cancer …/     residential history, NFHS/LASI working data
  TMC-BrCaProj/                 R analysis project (nested git repo)
    code/analysis/breast_cancer_risk_model.Rmd   the pipeline
    code/analysis/R/reference_imputation.R        hot-deck + reference builders
    data/reference/             generated reference panels
  presentations/                decks and running logs
  Papers/                       methodology literature
  docs/                         workspace documentation
```

---

## 16. Known limitations and caveats

- **Residence coding is 1 = urban, 2 = rural** across the model and reference
  (confirmed by the TMC, NFHS and LASI data dictionaries). Any consuming interface
  must send the same convention.
- **Single incidence table for both strata (by design).** Because the strata are
  age‑disjoint (Section 5), the total registry incidence maps to exactly one
  population at each age, so both models correctly share the one incidence table —
  this is a consequence of the age cutoff, not an approximation. A finer
  menopause‑specific λ_pre/λ_post split at the *same* age (via the age‑at‑menopause
  distribution) was prototyped and disabled; it would only be needed if the age
  bands were allowed to overlap. (Per‑individual scoring still applies total λ(t)
  at the woman’s actual age regardless of her self‑reported stratum.)
- **BMI × WHR interaction** is weakly supported (LRT p ≈ 0.11 pre, 0.06 post) and
  produces non‑monotone behaviour at extreme cells; BMI and WHR are collinear
  adiposity measures — candidate for simplification.
- **Imputed reference fields.** Breast lump and tobacco‑chewing timing are imputed
  (not observed) in NFHS/LASI.
- **BMI thresholds** are international (25, 30), not the lower Indian cutpoints.
- **Read‑but‑unused files.** Several 2023/2024 extracts are loaded but not merged
  (superseded by the Jan‑2026 drop); harmless but cluttering.
- **Mixed data vintages** in the model frame (Jan‑2026 anthro/PMH/chewing/smoking/
  alcohol; 2023 reproductive/demography/family history; Nov‑2025 residence; 2017
  questionnaire for menopause) — ID alignment confirmed via the funnel; the 2017
  and 2023 reproductive/menopause sources agree where they overlap (menopause
  exactly; parity/age‑at‑first‑pregnancy differ ~2–4.5%).
- **Phase‑2 / COVID waist‑to‑hip ratio.** Waist and hip circumference were not
  collected for Phase‑2 participants (enrolled during the COVID‑19 period, when
  contact anthropometry was suspended), so WHR is missing for ~36% of Phase‑2
  (≈21% of the analytic sample) and is multiply imputed for them; height and weight
  (BMI) were still available. The missingness is driven by enrolment period, not by
  WHR value or case status (missing‑at‑random conditional on phase), so hot‑deck
  imputation is defensible — but because WHR is a strong predictor and feeds the
  BMI×WHR interaction (caveat above), a complete‑WHR / Phase‑1 sensitivity analysis
  is recommended to confirm those effects are not imputation artifacts.

---

## 17. How to run

The pipeline is an R Markdown document. From the project root, with the required
packages installed (`iCARE`, `openxlsx`, `readstata13`, `haven`, `dplyr`,
`janitor`, `brglm2`, `naniar`, `purrr`):

```r
knitr::purl("TMC-BrCaProj/code/analysis/breast_cancer_risk_model.Rmd",
            output = "/tmp/pipeline.R")
source("/tmp/pipeline.R", echo = FALSE)
```

It reads the datasets in Section 2, builds the cohort (Section 3), fits the
models (Section 7), constructs the references (Section 8), computes absolute risk
(Section 10), runs the projection (Section 11) and validation (Section 12), and
saves the outputs in Section 13.
