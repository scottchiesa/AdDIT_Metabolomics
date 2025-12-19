************************************************************************************************************************************************
*******************************************************GRAPHING FIGURES FOR SUPPLEMENTAL FILE***************************************************
**************************************************************Author: Scott Chiesa**************************************************************
********************************************************************2023-2025*******************************************************************
************************************************************************************************************************************************

version 18
clear all
macro drop _all
set more off
set maxvar 10000

//T1D vs Controls//

cd "Y:\AdDIT\Metabolomics Analyses\Scott Analysis\Combined\"

import delimited "T1DvHC_supplementary_for_plot.csv", clear

rename v1 outcome
rename v2 beta
rename v3 lci
rename v4 uci
rename v5 p

gen row = _n

forvalues i = 1(41)205 {
    local j = `i' + 40
    local block = ceil(`i'/41)

    forestplot beta lci uci if inrange(row, `i', `j'), ///
        label(outcome) ///
        xline(0) xtitle("SD-difference in metabolite vs controls", size(vsmall)) ///
        title("") ///
		nostats ///
		xlabel(-1 1)

    graph export forest_block`block'.png, replace
}


//ACE vs Placebo//

import delimited "ACE_supplementary_for_plot.csv", clear

rename v1 outcome
rename v2 beta
rename v3 lci
rename v4 uci
rename v5 p

gen row = _n

forvalues i = 1(41)205 {
    local j = `i' + 40
    local block = ceil(`i'/41)

    forestplot beta lci uci if inrange(row, `i', `j'), ///
        label(outcome) ///
        xline(0) xtitle("SD-difference in metabolite vs placebo", size(vsmall)) ///
        title("") ///
		nostats ///
		xlabel(-1 1)

    graph export ACE_forest_block`block'.png, replace
}

//Statins vs Placebo//

import delimited "statin_supplementary_for_plot.csv", clear

rename v1 outcome
rename v2 beta
rename v3 lci
rename v4 uci
rename v5 p

gen row = _n

forvalues i = 1(41)205 {
    local j = `i' + 40
    local block = ceil(`i'/41)

    forestplot beta lci uci if inrange(row, `i', `j'), ///
        label(outcome) ///
        xline(0) xtitle("SD-difference in metabolite vs placebo", size(vsmall)) ///
        title("") ///
		nostats ///
		xlabel(-1 1)

    graph export statin_forest_block`block'.png, replace
}
