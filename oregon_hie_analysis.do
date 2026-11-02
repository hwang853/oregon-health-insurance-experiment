version 16.0
clear all
set more off

/*
Oregon Health Insurance Experiment analysis

Usage:
    stata -b do oregon_hie_analysis.do "/path/to/data" "results"

Both arguments are optional. The data directory defaults to the current
directory, and the output directory defaults to "results".
*/

args data_dir output_dir
if `"`data_dir'"' == "" local data_dir "."
if `"`output_dir'"' == "" local output_dir "results"

cd `"`data_dir'"'
capture mkdir `"`output_dir'"'

*******************************************************
***************        Problem 1        ***************
*******************************************************

* Problem 1 Question 1
use "oregonhie_descriptive_vars.dta", clear

merge 1:1 person_id using "oregonhie_stateprograms_vars.dta", assert(master match)
rename _merge _merge1
merge 1:1 person_id using "oregonhie_survey0m_vars.dta", assert(master match)
rename _merge _merge2
merge 1:1 person_id using "oregonhie_ed_vars.dta", assert(master match)
rename _merge _merge3

tab _merge3

* There are 74,922 observations in total.
* All individuals are observed in the descriptive, state-program, and
* baseline-survey files. Only 24,646 are observed in the ED file.

* Problem 1 Question 2
gen age = 2008 - birthyear_list
gen work_20 = (employ_hrs_0m == 3) | (employ_hrs_0m == 4)
replace work_20 = . if missing(employ_hrs_0m)

eststo: estpost summarize age female_list english_list work_20 er_any_0m, detail
esttab using `"`output_dir'/summary_statistics.tex"', replace ///
    cells("N mean(fmt(2)) sd(fmt(2)) p10th p50th p90th") noobs
eststo clear
* Problem 1 Question 3
gen treat = treatment

eststo: reg any_visit_ed treat numhh_list draw_lottery, cluster(household_id)
local ate = _b[treat]
display "ATE: `ate'"

esttab using `"`output_dir'/p1q3_ate.tex"', replace ///
    b(3) se(3) ar2(3) star(* 0.10 ** 0.05 *** 0.01) ///
    title("The ATE\label{tab1}") nonumbers mtitles("ER visit")

eststo clear

* Estimated ATE: 0.019. Medicaid increased the probability of an ER visit
* by about 1.9 percentage points; the estimate was significant at 1%.

* Problem 1 Question 4
count
local total_count = r(N)
count if treat == 1
local treat_count = r(N)

capture matrix drop ate_iterations
set seed 1234
forvalues i = 1/1000 {
    preserve
    quietly drop treat
    quietly gen rand = runiform()
    sort rand
    quietly gen treat = _n <= `treat_count'

    quietly reg any_visit_ed treat numhh_list draw_lottery
    matrix ate_iterations = nullmat(ate_iterations) \ _b[treat]
    restore
}

matrix list ate_iterations

* Temporarily replace the analysis data with the simulated estimates.
preserve
clear
svmat ate_iterations
summarize ate_iterations1, detail
count if ate_iterations1 >= `ate'
local fisher_pvalue = r(N) / 1000
display "Fisher P-value: `fisher_pvalue'"

histogram ate_iterations1, xline(`ate', lcolor(red)) ///
    graphregion(color(white)) ///
    xtitle("Placebo Treatment Effects of Oregon Health Program") ///
    color(gray%50)
graph export `"`output_dir'/random_inference_example.pdf"', replace
restore

*******************************************************
***************        Problem 2        ***************
*******************************************************

* Problem 2 Question 1
eststo: reg any_visit_ed ohp_all_ever_matchn_30sep2009 ///
    numhh_list draw_lottery, cluster(household_id)

esttab using `"`output_dir'/p2q1_result.tex"', replace ///
    b(3) se(3) ar2(3) star(* 0.10 ** 0.05 *** 0.01) ///
    title("The effect of Medicaid on ER visit\label{tab1}") ///
    nonumbers mtitles("ER visit")

eststo clear

* Problem 2 Question 2
gen agesqr = age^2

reg any_visit_ed ohp_all_ever_matchn_30sep2009 numhh_list draw_lottery ///
    age agesqr female_0m race_hisp_0m race_white_0m race_black_0m ///
    hhsize_0m edu_0m, cluster(household_id)

logit ohp_all_ever_matchn_30sep2009 numhh_list draw_lottery age agesqr ///
    female_0m race_hisp_0m race_white_0m race_black_0m hhsize_0m edu_0m

predict pscore

twoway (histogram pscore if ohp_all_ever_matchn_30sep2009 == 1, ///
        color(eltblue%40)) ///
    (histogram pscore if ohp_all_ever_matchn_30sep2009 == 0, ///
        color(erose%30)), ///
    legend(label(1 "Treated") label(2 "Control")) ///
    scheme(s1mono) xtitle("Propensity Score")

graph export `"`output_dir'/pscore_distribution.png"', replace

* Problem 2 Question 3
* Propensity-score matching
eststo ps_match: teffects psmatch (any_visit_ed) ///
    (ohp_all_ever_matchn_30sep2009 numhh_list draw_lottery age agesqr ///
     female_0m race_hisp_0m race_white_0m race_black_0m hhsize_0m edu_0m, ///
     logit), atet gen(ps_cps) nn(1)

* Inverse propensity-score weights
gen inv_ps_weight = ohp_all_ever_matchn_30sep2009 + ///
    (1 - ohp_all_ever_matchn_30sep2009) * pscore / (1 - pscore)

eststo ipw: reg any_visit_ed i.ohp_all_ever_matchn_30sep2009 ///
    [aw=inv_ps_weight], robust

esttab using `"`output_dir'/p2q3_result.tex"', replace ///
    b(3) se(3) ar2(3) star(* 0.10 ** 0.05 *** 0.01) ///
    title("The effect of Medicaid on ER visit\label{tab1}") ///
    nonumbers mtitles("PSM" "IPW")

eststo clear

