# Who Pays for a Tariff?

Border pass-through and rival-supplier pricing in the 2018–19 U.S.–China tariffs
using 2.45 million U.S. customs records (January 2017 – April 2019).

**Findings**
- U.S. importers paid the full tariff (pass-through 1.006); Chinese exporters did not cut prices.
- The 0.58 coefficient on statutory tariffs is a duty-collection rate, not partial pass-through.
- No reliable evidence that suppliers in other countries raised prices once China was tariffed.

## Data
FGKK replication package, Harvard Dataverse: https://doi.org/10.7910/DVN/KSOVSE
Download `m_flow_hs10_fm_new.dta`.

## Files
- `TARIFF_PROJECT_MASTER.R`: main replication and rival-supplier models (25 models).
- `referee_additions.R`: variation table, event study, lagged first stage, China-only check.
- `fig_event_study.png`: event study figure.

## Run
R 4.6.1 with fixest, data.table, haven, ggplot2. Set the data path at the top of each script.

    source("TARIFF_PROJECT_MASTER.R")
    results_folder <- run_project()

Author: Arthur Privette, UNC Chapel Hill · arthurprivette@gmail.com
