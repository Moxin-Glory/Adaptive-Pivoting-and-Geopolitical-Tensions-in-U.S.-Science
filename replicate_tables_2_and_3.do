************************************************************
* Replication do-file for Table 2, Table 3, and figure input data
* Data file: PNAS_revision_regression_upload_data.csv
* Outputs:
*   results/Table2_baseline_DID_10_columns.rtf
*   results/Table3_moderating_effects_6_columns.rtf
*   Fig3_event_study_data.csv
*   Fig4_did_topic_data.csv
************************************************************

version 17.0
clear all
set more off

************************************************************
* 0. Paths and packages
************************************************************

global input  "PNAS_revision_regression_upload_data.csv"
global outdir "results"
global event_csv "Fig3_event_study_data.csv"
global heterogeneity_csv "Fig4_did_topic_data.csv"

capture mkdir "$outdir"

foreach pkg in ftools reghdfe estout {
    capture which `pkg'
    if _rc ssc install `pkg', replace
}

import delimited "$input", clear varnames(1)

************************************************************
* 1. Check variables and define the panel
************************************************************

local required_vars ///
    author_id year treat ///
    log_us_federal_amount log_us_amount log_all_amount ///
    log_us_federal_grant_count log_us_grant_count ///
    log_all_grant_count reference_overlap ///
    word2vec_topic_exploration_30 concept_pivot_new_share_01 ///
    topic_change career_age log_paper_count log_cf ///
    log_all_us_collaborators_count ///
    ethnically_chinese high_list_related_field junior_scientist

foreach v of local required_vars {
    capture confirm variable `v'
    if _rc {
        display as error "Required variable not found: `v'"
        exit 111
    }
}

keep if inrange(year, 2014, 2023)

capture confirm numeric variable author_id
if _rc destring author_id, replace force

foreach v in author_fe post did did_topic sample_common ///
    es_2014 es_2015 es_2016 es_2018 es_2019 ///
    es_2020 es_2021 es_2022 es_2023 {
    capture drop `v'
}

capture isid author_id year
if _rc {
    display as error "author_id-year does not uniquely identify observations."
    duplicates report author_id year
    exit 459
}

egen long author_fe = group(author_id)
xtset author_fe year

************************************************************
* 2. DID variables, controls, and outcomes
************************************************************

gen byte post = (year >= 2019) if !missing(year)
gen byte did  = treat * post
gen double did_topic = did * topic_change

label variable did       "DID"
label variable did_topic "DID x Pivot Size"

global controls ///
    career_age ///
    log_paper_count ///
    log_cf ///
    log_all_us_collaborators_count

global y_fed          log_us_federal_amount
global y_us           log_us_amount
global y_total        log_all_amount
global y_pivot        topic_change
global y_fed_grants   log_us_federal_grant_count
global y_us_grants    log_us_grant_count
global y_total_grants log_all_grant_count
global y_word2vec     word2vec_topic_exploration_30
global y_concept      concept_pivot_new_share_01
global y_ref          reference_overlap

label variable career_age "Career Age"
label variable log_paper_count "# Papers"
label variable log_cf "Average Cf"
label variable log_all_us_collaborators_count "# U.S. Collaborators"
label variable topic_change "Pivot Size"

gen byte sample_common = !missing( ///
    author_fe, year, treat, did, ///
    career_age, log_paper_count, log_cf, ///
    log_all_us_collaborators_count)

************************************************************
* 3. Table 2: baseline DID regressions
************************************************************

eststo clear

local table2_y1  "$y_fed"
local table2_y2  "$y_us"
local table2_y3  "$y_total"
local table2_y4  "$y_pivot"
local table2_y5  "$y_fed_grants"
local table2_y6  "$y_us_grants"
local table2_y7  "$y_total_grants"
local table2_y8  "$y_word2vec"
local table2_y9  "$y_concept"
local table2_y10 "$y_ref"

forvalues column = 1/10 {
    local yvar "`table2_y`column''"

    reghdfe `yvar' did $controls i.year ///
        if sample_common == 1 & !missing(`yvar'), ///
        absorb(author_fe) vce(cluster author_fe)

    eststo table2_`column'
    estadd local YearFE "Y", replace
    estadd local IndividualFE "Y", replace
}

esttab table2_1 table2_2 table2_3 table2_4 table2_5 ///
    table2_6 table2_7 table2_8 table2_9 table2_10 ///
    using "$outdir/Table2_baseline_DID_10_columns.rtf", ///
    replace nogaps nonotes b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.001) ///
    keep(did career_age log_paper_count log_cf ///
        log_all_us_collaborators_count _cons) ///
    order(did career_age log_paper_count log_cf ///
        log_all_us_collaborators_count _cons) ///
    coeflabels( ///
        did "DID" ///
        career_age "Career Age" ///
        log_paper_count "# Papers" ///
        log_cf "Average Cf" ///
        log_all_us_collaborators_count "# U.S. Collaborators" ///
        _cons "Constant") ///
    mtitles( ///
        "U.S. Federal Funding" ///
        "U.S. Funding" ///
        "Total Funding" ///
        "Pivot Size" ///
        "# U.S. Federal Grants" ///
        "# U.S. Grants" ///
        "# Total Grants" ///
        "Word2Vec Pivot Size" ///
        "Concept-based Pivot Size" ///
        "Reference Overlap") ///
    stats(YearFE IndividualFE N r2_a, ///
        labels("Year FE" "Individual FE" "N" "Adj. R2") ///
        fmt(%s %s %12.0fc %9.3f)) ///
    addnotes( ///
        "Sample sizes vary across models because non-missing outcome values are required." ///
        "Standard errors, clustered at the individual level, are reported in parentheses." ///
        "*** p < 0.001, ** p < 0.05, * p < 0.1.") ///
    title("Table 2. The impact of U.S.-China geopolitical tensions on funding amount and pivot size.")

************************************************************
* 4. Table 3: pivot-size moderation regressions
************************************************************

eststo clear

local table3_y1 "$y_fed"
local table3_y2 "$y_us"
local table3_y3 "$y_total"
local table3_y4 "$y_fed_grants"
local table3_y5 "$y_us_grants"
local table3_y6 "$y_total_grants"

forvalues column = 1/6 {
    local yvar "`table3_y`column''"

    reghdfe `yvar' did_topic did topic_change $controls i.year ///
        if sample_common == 1 & !missing(`yvar', topic_change), ///
        absorb(author_fe) vce(cluster author_fe)

    eststo table3_`column'
    estadd local YearFE "Yes", replace
    estadd local IndividualFE "Yes", replace
}

esttab table3_1 table3_2 table3_3 table3_4 table3_5 table3_6 ///
    using "$outdir/Table3_moderating_effects_6_columns.rtf", ///
    replace nogaps nonotes b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.001) ///
    keep(did_topic did topic_change career_age log_paper_count ///
        log_cf log_all_us_collaborators_count _cons) ///
    order(did_topic did topic_change career_age log_paper_count ///
        log_cf log_all_us_collaborators_count _cons) ///
    coeflabels( ///
        did_topic "DID x Pivot Size" ///
        did "DID" ///
        topic_change "Pivot Size" ///
        career_age "Career Age" ///
        log_paper_count "# Papers" ///
        log_cf "Average Cf" ///
        log_all_us_collaborators_count "# U.S. Collaborators" ///
        _cons "Constant") ///
    mtitles( ///
        "U.S. Federal Funding" ///
        "U.S. Funding" ///
        "Total Funding" ///
        "# U.S. Federal Grants" ///
        "# U.S. Grants" ///
        "# Grants") ///
    stats(YearFE IndividualFE N r2_a, ///
        labels("Year FE" "Individual FE" "N" "Adj. R2") ///
        fmt(%s %s %12.0fc %9.3f)) ///
    addnotes( ///
        "*** p < 0.001, ** p < 0.05, * p < 0.1." ///
        "Standard errors, clustered at the individual level, are reported in parentheses.") ///
title("Table 3. Associations between pivot size and differential changes in funding.")

************************************************************
* 5. Figure 2 input data: event-study estimates
************************************************************

foreach y in 2014 2015 2016 2018 2019 2020 2021 2022 2023 {
    gen byte es_`y' = (year == `y' & treat == 1)
    label variable es_`y' "`y' x Treat"
}

tempfile event_results
tempname event_handle

postfile `event_handle' str32 outcome int panel int year ///
    double estimate se ci_low ci_high using `event_results', replace

local panelno = 0
foreach outcome in fed us total pivot {
    local ++panelno
    if "`outcome'" == "fed" {
        local yvar "$y_fed"
        local outname "Federal Funding"
    }
    if "`outcome'" == "us" {
        local yvar "$y_us"
        local outname "U.S. Funding"
    }
    if "`outcome'" == "total" {
        local yvar "$y_total"
        local outname "Total Funding"
    }
    if "`outcome'" == "pivot" {
        local yvar "$y_pivot"
        local outname "Pivot Size"
    }

    reghdfe `yvar' ///
        es_2014 es_2015 es_2016 es_2018 es_2019 ///
        es_2020 es_2021 es_2022 es_2023 ///
        $controls i.year ///
        if sample_common == 1 & !missing(`yvar'), ///
        absorb(author_fe) vce(cluster author_fe)

    scalar crit = invttail(e(df_r), 0.025)

    forvalues y = 2014/2023 {
        if `y' == 2017 {
            post `event_handle' ("`outname'") (`panelno') (`y') (0) (0) (0) (0)
        }
        else {
            capture scalar b_event = _b[es_`y']
            if _rc {
                scalar b_event = .
                scalar se_event = .
                scalar lo_event = .
                scalar hi_event = .
            }
            else {
                scalar se_event = _se[es_`y']
                scalar lo_event = b_event - crit * se_event
                scalar hi_event = b_event + crit * se_event
            }
            post `event_handle' ("`outname'") (`panelno') (`y') ///
                (b_event) (se_event) (lo_event) (hi_event)
        }
    }
}

postclose `event_handle'

preserve
use `event_results', clear
sort panel year
assert _N == 40
assert estimate == 0 & ci_low == 0 & ci_high == 0 if year == 2017
export delimited using "$event_csv", replace
restore

************************************************************
* 6. Figure 3 input data: heterogeneity estimates
************************************************************

foreach g in ethnically_chinese high_list_related_field junior_scientist {
    assert inlist(`g', 0, 1, .)
}

eststo clear
local model = 0

foreach outcome in fed us {
    if "`outcome'" == "fed" local yvar "$y_fed"
    if "`outcome'" == "us"  local yvar "$y_us"

    foreach groupvar in ethnically_chinese high_list_related_field junior_scientist {
        foreach groupvalue in 1 0 {
            local ++model

            reghdfe `yvar' did_topic did topic_change $controls i.year ///
                if sample_common == 1 ///
                & `groupvar' == `groupvalue' ///
                & !missing(`yvar', topic_change), ///
                absorb(author_fe) vce(cluster author_fe)

            eststo m`model'
        }
    }
}

assert `model' == 12

esttab m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 ///
    using "$heterogeneity_csv", ///
    replace keep(did_topic) ///
    cells("b(fmt(3)) p(fmt(3)) ci_l(fmt(3)) ci_u(fmt(3))") ///
    compress nonumber csv ///
    mtitles( ///
        "Federal: Ethnically Chinese scientists" ///
        "Federal: Non-ethnically Chinese scientists" ///
        "Federal: High list-related fields" ///
        "Federal: Low list-related fields" ///
        "Federal: Junior scientists" ///
        "Federal: Senior scientists" ///
        "U.S.: Ethnically Chinese scientists" ///
        "U.S.: Non-ethnically Chinese scientists" ///
        "U.S.: High list-related fields" ///
        "U.S.: Low list-related fields" ///
        "U.S.: Junior scientists" ///
        "U.S.: Senior scientists")

display as result "Done. Table outputs were saved in: $outdir"
display as result "Done. Figure input data were saved as: $event_csv and $heterogeneity_csv"
