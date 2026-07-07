# NHANES VOC Demo: DAG, Variables, and VOC Groups

This note summarizes only the causal structure, variable list, and exposure
grouping for the NHANES VOC method demonstration.

## DAG

Target contrast:

- Exposure: urinary VOC metabolite mixture
- Outcome: log-transformed gamma-glutamyl transferase, `log_ggt`
- Interpretation: higher `log_ggt` is less favorable; the prespecified primary
  direction is positive.

DAGitty specification:

```r
dag {
Age [pos="0,2"]
Sex [pos="0,1"]
RaceEthnicity [pos="0,0"]
SES [pos="1,0"]
Smoking [pos="1,2"]
BMI [pos="2,2"]
UrineDilution [pos="2,0"]
VOC [exposure,pos="3,1"]
GGT [outcome,pos="5,1"]

Age -> BMI
Age -> GGT
Age -> UrineDilution
Age -> VOC
Sex -> BMI
Sex -> GGT
Sex -> UrineDilution
Sex -> VOC
RaceEthnicity -> GGT
RaceEthnicity -> SES
RaceEthnicity -> Smoking
RaceEthnicity -> VOC
SES -> GGT
SES -> Smoking
SES -> VOC
Smoking -> GGT
Smoking -> VOC
BMI -> GGT
BMI -> UrineDilution
BMI -> VOC
UrineDilution -> VOC
VOC -> GGT
}
```

DAG-based adjustment variables:

- Age
- Sex
- Race/ethnicity
- Socioeconomic status
- Smoking
- BMI

Urinary creatinine is added separately as a dilution adjustment rather than as
a backdoor confounder.

## Variable List

Outcome:

| Role | NHANES variable | Model variable | Notes |
|---|---|---|---|
| Outcome | `LBXSGTSI` | `log_ggt` | `log(LBXSGTSI)`; higher values are less favorable |

Survey design:

| Role | NHANES variable | Notes |
|---|---|---|
| PSU | `SDMVPSU` | `ids = ~SDMVPSU` |
| Strata | `SDMVSTRA` | `strata = ~SDMVSTRA` |
| VOC subsample weight | `WTSAPRP` | `weights = ~WTSAPRP`, for P_UVOC |

DAG and model covariates:

| DAG node | NHANES variable | Model variable | Notes |
|---|---|---|---|
| Age | `RIDAGEYR` | `age` | Continuous age |
| Sex | `RIAGENDR` | `sex` | Male/Female factor |
| Race/ethnicity | `RIDRETH3` | `race` | Main race/ethnicity factor |
| Race/ethnicity | `RIDRETH3` | `race4` | Collapsed sensitivity variable |
| SES | `INDFMPIR` | `pir` | Poverty income ratio |
| Smoking | `LBXCOT` | `log_cotinine` | `log1p(LBXCOT)` |
| BMI | `BMXBMI` | `bmi` | Body mass index |
| Urine dilution | `URXUCR` | `log_ucr` | `log(URXUCR)` |

Implemented covariate sets:

| Set | Variables |
|---|---|
| `minimal` | `age`, `sex`, `race` |
| `ucr` | `age`, `sex`, `race`, `log_ucr` |
| `bmi_ucr` | `age`, `sex`, `race`, `bmi`, `log_ucr` |
| `dag_primary` | `age`, `sex`, `race`, `pir`, `log_cotinine`, `bmi`, `log_ucr` |
| `dag_no_tobacco` | `age`, `sex`, `race`, `pir`, `bmi`, `log_ucr` |
| `survey_reduced` | `age`, `sex`, `race4`, `bmi`, `log_ucr` |

## VOC Groups

The exposure set uses 16 P_UVOC analytes with unweighted detection rate at or
above 75%; `URXBPM` is excluded.

| Group | Variable | Metabolite | Parent/proxy exposure |
|---|---|---|---|
| Aromatic hydrocarbons | `URX2MH` | 2-methylhippuric acid | o-xylene |
| Aromatic hydrocarbons | `URX34M` | 3-/4-methylhippuric acid | m-/p-xylene |
| Aromatic hydrocarbons | `URXBMA` | N-acetyl-S-(benzyl)-L-cysteine | toluene |
| Aromatic hydrocarbons | `URXMAD` | Mandelic acid | styrene/ethylbenzene |
| Aromatic hydrocarbons | `URXPHG` | Phenylglyoxylic acid | styrene/ethylbenzene |
| Dienes/alkenes/epoxides | `URXDHB` | N-acetyl-S-(3,4-dihydroxybutyl)-L-cysteine | 1,3-butadiene |
| Dienes/alkenes/epoxides | `URXMB3` | N-acetyl-S-(4-hydroxy-2-butenyl)-L-cysteine | 1,3-butadiene |
| Dienes/alkenes/epoxides | `URXIPM3` | N-acetyl-S-(4-hydroxy-2-methyl-2-buten-1-yl)-L-cysteine | isoprene |
| Dienes/alkenes/epoxides | `URXHP2` | N-acetyl-S-(2-hydroxypropyl)-L-cysteine | propylene oxide |
| Nitrile/amide/cyanide-related | `URXAAM` | N-acetyl-S-(2-carbamoylethyl)-L-cysteine | acrylamide |
| Nitrile/amide/cyanide-related | `URXAMC` | N-acetyl-S-(N-methylcarbamoyl)-L-cysteine | N,N-dimethylformamide |
| Nitrile/amide/cyanide-related | `URXCEM` | N-acetyl-S-(2-carboxyethyl)-L-cysteine | acrylonitrile/acrolein-related |
| Nitrile/amide/cyanide-related | `URXCYM` | N-acetyl-S-(2-cyanoethyl)-L-cysteine | acrylonitrile |
| Nitrile/amide/cyanide-related | `URXATC` | 2-aminothiazoline-4-carboxylic acid | cyanide-related exposure |
| Reactive aldehydes | `URXHPM` | N-acetyl-S-(3-hydroxypropyl)-L-cysteine | acrolein |
| Reactive aldehydes | `URXPMM` | N-acetyl-S-(3-hydroxypropyl-1-methyl)-L-cysteine | crotonaldehyde |

