# VaccineValidateR

**VaccineValidateR** is an R package for testing vaccine histories against the ACIP immunization recommendations.
It is designed for **research and data quality purposes** — not for clinical decision-making.  

**Disclaimer:** This package is *not* an official ACIP/CDC tool and should not be used for clinical decision support.  
Use is limited to research, evaluation, and educational contexts.

---

## Features

- Converts immunization data into antigen-based vaccination series (using a CVX code to antigen map).
  - CVX codes updated to 09/01/2025
  - Handles some non-US immunization products without CVX codes
- Vaccination validation
  - Applies age and interval rules from most recent valid [ACIP](https://www.cdc.gov/vaccines/hcp/imz-schedules/index.html) (updated 09/01/2025) guidance for each antigen/dose.
  - Considers edge cases and situational handling as per guidance on [Immunize.org](https://www.immunize.org/)
  - Addresses numerous special cases (examples):
    - Live vaccine minimum spacing across multiple series
    - Meningitis B product-specific timing calculations
  - Evaluates delayed delivery dose-by-dose (e.g. if a first dose is delayed, calculation of any subsequent dose delays uses the date of the first administration to assess for delay)
  - Flags series completion for main childhood antigens
  - Provides three main outputs:
    - **immunizations** (product-level data mapping back content to input)
    - **antigens** (dose-level validated data)
    - **skipped_antigens** (administrations that do not meet ACIP timing - some of these are interval doses that are ignored or doses given for travel or as part of non-US schedules)
- Visit level evaluations
  - *Used after validation* 
  - Identifies missed opportunities and whether immunizations given at or after a visit were delayed
  - Notes if a future dose was given wthin 15, 30, or 90 days after a missed opportunity
  - Considers each combination of visit and antigen that could be given at a visit
- Summarization functions 
  - * Used after validation*
  - Population level
    - Computes composite completion measures (e.g., **HEDIS** child/adolescent metrics, UTD status evaluations, special popluation metrics).
  - Visit level
    - Aggregates visit-level results by antigen (counts, % missed, % caught up within X days), grouped by year, antigen, or custom grouping columns (e.g. SYSTEM/VISIT_TYPE).
- Extensible
  - Designed to support new ACIP changes, additional antigens, and site/network-level QI.
  - Scales to large datasets (tested on >60M administrations)
---

## Installation

This package is not on CRAN.  
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
```
### Status

This package is updated in response to changes in ACIP rules yearly (and sometimes more frequently).
Bug/error fixes are released when noticed.
Validation rules are based on ACIP 2025 guidance, but may not yet cover every nuance.
COVID vaccination status is not calculated for completion given the recent changes to yearly formulations.
New features to be added as time permits. Ideas welcome.

### Acknowledgements
The author wished to acknowledge colleages and early users for their feedback, guidance, and suggestions including Dr. Alex Fiks, Dr. Robert Grundmeier, Mary Kate Kelly, and Abbie Steiner.

### Non-Commercial Research License
---
This project is released under a Non-Commercial Research License. For commercial use, please contact michelj@chop.edu for licensing terms.
Copyright ©2025 The Children's Hospital of Philadelphia. See full licenses statement for more details.
Permission is hereby granted, free of charge, to any person or organization to use, copy, modify, and distribute this software and associated documentation files (the “Software”), for academic, research, or educational purposes only, subject to the following conditions:

#### Attribution
Appropriate credit must be given to the authors in any use, publication, or derivative work of the Software.

#### Non-Commercial Use Only
The Software may not be used, in whole or in part, for commercial purposes, including but not limited to:
- use in a product for sale,
- use in a for-profit company’s operations,
- use in services provided to customers for a fee.

#### Commercial Licensing
For commercial use, a separate license must be obtained from the copyright holder. Please contact:
[michelj@chop.edu](mailto:michelj@chop.edu)

#### Warranty Disclaimer
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Resources

- [ACIP Immunization Schedules (CDC)](https://www.cdc.gov/vaccines/hcp/imz-schedules/index.html)
- [CDC CVX Code List](https://www.cdc.gov/vaccines/programs/iis/codes.html)
- [devtools GitHub](https://github.com/r-lib/devtools)

### Dose Calculation Tables (For Reference)

These tables summarize **minimum ages, intervals, and series-completion rules** for common antigens.
Detailed dose timing and series completion tables are provided in the [dose reference vignette](/docs/articles/dose_reference.html).
