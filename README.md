# VaccineValidateR

**VaccineValidateR** is an R package for testing vaccine histories against the ACIP immunization recommendations.
It is designed for **research and data quality purposes** — not for clinical decision-making.  

**Disclaimer:** This package is *not* an official ACIP/CDC tool and should not be used for clinical decision support.  
Use is limited to research, evaluation, and educational contexts.

---

## Features

- Preprocesses raw immunization data (CVX codes → antigens).
- Applies age and interval rules from [ACIP](https://www.cdc.gov/vaccines/hcp/imz-schedules/index.html)(updated 09/01/2025) for each antigen/dose.
- Edge cases and specific situation handling as per guidance on [Immunize.org](https://www.immunize.org/)
- CVX codes updated to 09/01/2025
- Addresses numerous special cases (examples):
  - Live vaccine minimum spacing between series
  - Meningitis B product selection
- Scales to very large datasets (tested on >60M administrations)
- Handles some non-US immunization products without CVX codes
- Validation function:
  - Flags series completion
  - Evaluates delayed delivery dose-by-dose (e.g. if a first dose is delayed, calculation of any subsequent dose delays uses the date of the first administration to assess for delay)
  - Provides four main outputs:
    - **immunizations** (product-level data)
    - **antigens** (dose-level validated data)
    - **invalid_immunizations** (unmapped or unprocessable records)
    - **invalid_antigens** (administrations that do not meet ACIP timing - some of these are interval doses that are ignored or doses given for travel or as part of non-US schedules)
- Evaluate_visits
  - *Used after validation* 
  - Determines missed opportunities
  - Calculates if future doses were delayed after missed opportunities
- Population_metrics function 
  - * Used after validation*
  - Calculates completion of multiple HEDIS measures

---

## Installation

This package is not yet on CRAN.  
You can install directly from GitHub using:

```r
# install.packages("devtools")
devtools::install_github("psyllix/VaccineValidateR")

### Example
r

library(VaccineValidateR)

---

# Example immunization data
data <- data.table::data.table(
  STUDY_ID = 1,
  CVX = c("8", "8"),#HEPATITIS B
  AGE_IMM_GIVEN = c("0","31"),
  DATE_GIVEN = as.Date(c("2022-01-01", "2022-02-01"))#BIRTH AND 1 MONTH DOSE
)

# Validate
result <- VaccineValidateR::validation(data)

# Check antigen-level results
head(result$antigens)

### Status

This package is updated in response to changes in ACIP rules yearly (and sometimes more frequently).
Bug/error fixes are released when notices.
Validation rules are based on ACIP 2024 guidance, but may not yet cover every nuance.
COVID vaccination status is not calculated for completion given the recent changes to yearly formulations.

### License
---

This project is licensed under the [GNU General Public License v3 (GPL-3)](https://www.gnu.org/licenses/gpl-3.0.html).
You are free to use, modify, and distribute this software, provided that derivative works are also licensed under GPL-3.
See the GNU General Public License for details.

## Resources

- [ACIP Immunization Schedules (CDC)](https://www.cdc.gov/vaccines/hcp/imz-schedules/index.html)
- [CDC CVX Code List](https://www.cdc.gov/vaccines/programs/iis/codes.html)
- [devtools GitHub](https://github.com/r-lib/devtools)

### Dose Calculation Tables (For Reference)

| Antigen    | Dose Number | Minimum Age | Minimum Interval | Notes                                                                 | Completes Series |
|------------|------------|------------|----------------|----------------------------------------------------------------------|----------------|
| COVID      | 1          | 6 months   | NA             | None                                                                 | No             |
| HEPA       | 1          | 1 year     | NA             | None                                                                 | No             |
| HEPB       | 1          | 0 years    | NA             | None                                                                 | No             |
| HIB        | 1          | 15 months  | NA             | None                                                                 | Yes            |
| HIB        | 1          | 6 weeks    | NA             | None                                                                 | No             |
| HPV        | 1          | 9 years    | NA             | None                                                                 | No             |
| INFLUENZA  | 1          | 6 months   | NA             | None                                                                 | No             |
| MCV        | 1          | 16 years   | NA             | None                                                                 | Yes            |
| MCV        | 1          | 10 years   | NA             | None                                                                 | No             |
| MENB       | 1          | 16 years   | NA             | Start of SCDM series cannot be delayed                               | No             |
| MMR        | 1          | 1 year     | NA             | None                                                                 | No             |
| PCV        | 1          | 2 years    | NA             | None                                                                 | Yes            |
| PCV        | 1          | 6 weeks    | NA             | None                                                                 | No             |
| POLIO      | 1          | 6 weeks    | NA             | None                                                                 | No             |
| ROTA       | 1          | 4 weeks    | NA             | Max age to start: 14w6d                                              | No             |
| RSV        | 1          | 0 years    | NA             | None                                                                 | Yes            |
| TETANUS    | 1          | 6 weeks    | NA             | None                                                                 | No             |
| VZV        | 1          | 1 year     | NA             | None                                                                 | No             |

| Antigen    | Dose Number | Minimum Age | Minimum Interval | Notes                                                                 | Completes Series |
|------------|------------|------------|----------------|----------------------------------------------------------------------|----------------|
| COVID      | 2          | 6 months   | 3–4 wks        | Interval depends on first dose brand (Pfizer 3 wk, others 4 wk)      | No             |
| HEPA       | 2          | 1–19 years | 6 months       | None                                                                 | Yes            |
| HEPA       | 2          | Adult >19 yr | 6 months     | None                                                                 | Yes            |
| HEPA       | 2          | Adult >19 yr | 4 wks        | None                                                                 | No             |
| HIB        | 2          | 15 months  | 8 wks          | None                                                                 | Yes            |
| HIB        | 2          | 1 year     | 8 wks          | First dose ≥1 yr                                                     | Yes            |
| HIB        | 2          | 10 weeks   | 4 wks          | None                                                                 | No             |
| HPV        | 2          | 9 years    | 5 months       | First dose <15 yr                                                    | Yes            |
| HPV        | 2          | 9 years    | 4 wks          | None                                                                 | No             |
| INFLUENZA  | 2          | 6 months   | 1 year         | Yearly booster; cannot be delayed                                     | No             |
| MENB       | 2          | 16 years   | 1 month        | Both doses prior to date of change in guidance | Yes            |
| MENB       | 2          | 16 years   | 6 months       | Both doses Bexsero                                                    | Yes            |
| MENB       | 2          | 16 years   | 4 months       | Both doses Bexsero* | No       |
| MENB       | 2          | 16 years   | 6 months       | Both doses Trumenba                                                   | Yes            |
| MMR        | 2          | 1 year     | 4 wks          | None                                                                 | Yes            |
| PCV        | 2          | 1 year     | 8 wks          | None                                                                 | No             |
| PCV        | 2          | 2 years    | 8 wks          | None                                                                 | Yes            |
| PCV        | 2          | 10 wk–1 yr | 4 wks          | None                                                                 | No             |
| POLIO      | 2          | 10 weeks   | 4 wks          | None                                                                 | No             |
| ROTA       | 2          | 10 weeks   | 4 wks          | Both first and second dose are 2-dose rotavirus                       | Yes            |
| ROTA       | 2          | 10 weeks   | 4 wks          | None                                                                 | No             |
| TETANUS    | 2          | 10 weeks   | 4 wks          | None                                                                 | No             |
| VZV        | 2          | 1 year     | 12 wks         | None                                                                 | Yes            |
| VZV        | 2          | 13 years   | 4 wks          | None                                                                 | Yes            |

| Antigen    | Dose Number | Minimum Age | Minimum Interval | Notes                                                                 | Completes Series |
|------------|------------|------------|----------------|----------------------------------------------------------------------|----------------|
| COVID      | 3          | 6 months   | 4–8 wks        | Interval depends on first dose brand: Pfizer 8 wks, others 4 wks      | Yes            |
| HEPA       | 3          | 19+ years  | 5 months       | Interval from first dose ≥6 mo                                         | Yes            |
| HEPB       | 3          | 24 weeks   | 8 wks          | 16 wks from dose 1                                                    | Yes            |
| HIB        | 3          | 1 year     | 8 wks          | All previous doses are HIB3                                           | Yes            |
| HIB        | 3          | 12–24 months | 8 wks        | First dose ≥7 mo                                                      | Yes            |
| HIB        | 3          | 12+ months | 8 wks          | Last dose ≥12 mo                                                      | Yes            |
| HIB        | 3          | 14 weeks   | 4 wks          | None                                                                 | No             |
| HPV        | 3          | 9 years    | 12 wks         | Interval from first dose ≥21 wks                                       | Yes            |
| PCV        | 3          | 24 months  | 8 wks          | None                                                                 | Yes            |
| PCV        | 3          | 12+ months | 8 wks          | None                                                                 | Yes            |
| PCV        | 3          | 14 weeks   | 4 wks          | None                                                                 | No             |
| POLIO      | 3          | 4 years    | 180 days       | None                                                                 | Yes            |
| POLIO      | 3          | 14 weeks   | 4 wks          | None                                                                 | No             |
| ROTA       | 3          | 14 weeks   | 4 wks          | None                                                                 | Yes            |
| TETANUS    | 3          | 10 years   | 180 days       | Adolescent Vaccination                                                | Never Completed|
| TETANUS    | 3          | 14 weeks   | 4 wks          | None                                                                 | No             |

| Antigen    | Dose Number | Minimum Age | Minimum Interval | Notes                                                                 | Completes Series |
|------------|------------|------------|----------------|----------------------------------------------------------------------|----------------|
| HIB        | 4          | 1 year     | 8 wks          | None                                                                 | Yes            |
| PCV        | 4          | 1 year     | 8 wks          | None                                                                 | Yes            |
| POLIO      | 4          | 4 years    | 180 days       | None                                                                 | Yes            |
| TETANUS    | 4          | 10 years   | 180 days       | Adolescent Vaccination                                                | Never Completed|
| TETANUS    | 4          | 4–6 years  | 180 days       | None                                                                 | No             |
| TETANUS    | 4          | 7+ years   | 180 days       | None                                                                 | No             |
| TETANUS    | 4          | 1–6 years  | 4 months       | None                                                                 | No             |

| Antigen    | Dose Number | Minimum Age | Minimum Interval | Notes                                                                 | Completes Series |
|------------|------------|------------|----------------|----------------------------------------------------------------------|----------------|
| TETANUS    | 5          | 10 years   | 180 days       | Adolescent Vaccination                                                | Never Completed|
| TETANUS    | 5          | 4–6 years  | 180 days       | None                                                                 | No             |
| TETANUS    | 6          | 10 years   | 180 days       | Adolescent Vaccination                                                | Never Completed|
| TETANUS    | 7          | 10 years   | –              | Adolescent Vaccination                                                | Never Completed|

| Antigen    | Dose Number | Minimum Age | Minimum Interval | Notes                                                                 | Completes Series |
|------------|------------|------------|----------------|----------------------------------------------------------------------|----------------|
| COVID      | Yearly     | 6 months   | 8 weeks        | 2024+ logic                                                           | No             |
| INFLUENZA  | Yearly     | –          | –              | Yearly dose                                                           | No             |

*
