*************************************************************************************************************************************************
****************************************************************SUPPLEMENTARY ANALYSIS***********************************************************
*****************************************************************Author: Scott Chiesa************************************************************
***********************************************************************2023-2025*****************************************************************
*************************************************************************************************************************************************

version 18
clear all
macro drop _all
set more off
set maxvar 10000

cd "Y:\AdDIT\Metabolomics Analyses\Scott Analysis\Combined"
use "Master Dataset.dta", clear

set matsize 11000

//Drop some weird rows that have exported over in dataset//

drop in 785/1042

//Fix missing values in dataset//
 
replace Country = 2 in 704/725
replace Country = 3 in 726/784
 
//Rename variables//

rename sex2 sex
rename ageatvisit age
rename durationatvisit duration
rename Height height
rename Weight weight
rename Waist waist
rename BMI bmi
rename SBPmean sbp
rename DBPmean dbp
rename HbA1cperc hba1c
rename H_L_Risk h_l_risk
rename ACE_v_Placebo a_v_p
rename Placebo_v_Obs p_v_o
rename Statin_v_Placebo s_v_p
rename HDLChol hdl
rename LDLChol ldl
rename Trigl trig
rename hsCRP crp
rename Country country

//Generate Baseline T1D vs Control Comparison Groups//
 
gen DvC = .
replace DvC = 0 if visit == 2
replace DvC = 1 if visit == 0

//Put variables in order to allow loops to be run for all outcomes from here on//

order visit barcode sex age ageatdiagnosis duration height weight bmi waist sbp dbp ldl hdl trig hba1c country h_l_risk a_v_p s_v_p p_v_o 
drop Total_C Remnt_C Total_TG Total_CE Total_FC VLDL_C VLDL_TG LDL_C LDL_TG HDL_C HDL_TG ApoB ApoA1 ApoB_by_ApoA1 Phosphoglyc Cholines Phosphatidylc Sphingomyelins Total_FA SFA MUFA DHA LA PUFA Omega_3 Omega_6 Ala Gln His Phe Tyr Ile Leu Val Glucose Lactate Pyruvate Citrate Glycerol Acetate Acetoacetate bOHbutyrate Creatinine Albumin GlycA

//Log transform everything for later analyses looking at % change//

foreach var of varlist(non_HDL_C - S_HDL_TG_pct) {
	gen log`var' = log(`var')
}

//z-score all potential metabolomic outcomes for comparison//

zscore non_HDL_C - S_HDL_TG_pct

//Run baseline comparisons of people with vs without diabetes//

 ds z_non_HDL_C - z_S_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolite: `metabolite'"
  eststo: quietly mixed `metabolite' DvC age sex bmi || country:
  }
  esttab _all using T1DvHC_supplementary.csv, replace cells(b(fmt(3)) se(fmt(3)) p(fmt(3))) gaps lines nostar keep(DvC)
eststo clear

//Put in nice order//

egen patient = group(barcode)
xtset patient visit
order patient visit
sort patient visit


****************************************CREATING BASELINE AND FINAL VARIBALES FOR OUTCOMES AND ADJUSTMENTS***************************************


//Metabolites//

ds non_HDL_C - S_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolite: `metabolite'"
  by patient (visit), sort: gen bl`metabolite' = `metabolite'[1]
by patient (visit), sort: gen final`metabolite' = `metabolite'[2]
  }
  
//Log Metabolites//

ds lognon_HDL_C - logS_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolite: `metabolite'"
  by patient (visit), sort: gen bl`metabolite' = `metabolite'[1]
by patient (visit), sort: gen final`metabolite' = `metabolite'[2]
  }

//Metabolite z-scores//

ds z_non_HDL_C - z_S_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolite: `metabolite'"
  by patient (visit), sort: gen bl`metabolite' = `metabolite'[1]
by patient (visit), sort: gen final`metabolite' = `metabolite'[2]
  }

//Other covariates//  
  
by patient (visit), sort: gen blbmi = bmi[1]
by patient (visit), sort: gen finalbmi = bmi[2]
by patient (visit), sort: gen blsbp = sbp[1]
by patient (visit), sort: gen bldbp = dbp[1]
by patient (visit), sort: gen blage = age[1]
by patient (visit), sort: gen finalage = age[2]
by patient (visit), sort: gen blduration = duration[1]
by patient (visit), sort: gen finalduration = duration[2]
by patient (visit), sort: gen blhba1c = hba1c[1]
by patient (visit), sort: gen finalhba1c = hba1c[2]

//create time variable to account for different times on drugs etc//

gen time = finalage-blage 

//create indicator for those followed observationally//

gen notreat = 0
replace notreat = 1 if treatgroups == 4
replace notreat = 1 if treatgroups == 5


***************************************************COMPARISON OF DRUG EFFECTS IN RCT PATIENTS ONLY***********************************************
********************************(*********************REGRESSION MODELS COMPARING FINAL VALUES***************************************************

//ACE v Placebo//

recode a_v_p (1=4) (2=3)
recode a_v_p (3=0) (4=1)

  ds z_non_HDL_C - z_S_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolites: `metabolite'"
  eststo: quietly regress `metabolite' a_v_p blage sex blduration blhba1c bl`metabolite' country if visit ==1,  
  }
  esttab _all using ACEfinalresults_streamlined_supplementary.csv, replace cells(b(fmt(3)) se(fmt(3)) p(fmt(3))) gaps lines nostar keep(a_v_p)
  eststo clear
  
//Statin v Placebo//

recode s_v_p (1=4) (2=3)
recode s_v_p (3=0) (4=1)

  ds z_non_HDL_C - z_S_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolites: `metabolite'"
  eststo: quietly regress `metabolite' s_v_p blage sex blduration blhba1c bl`metabolite' country if visit ==1,  
  }
  esttab _all using statinfinalresults_streamlined_supplementary.csv, replace cells(b(fmt(3)) se(fmt(3)) p(fmt(3))) gaps lines nostar keep(s_v_p)
  eststo clear
  

***************************************************LONGITUDINAL CHANGE IN T1D OVER TIME***************************************************

  
//Drop any treated patients as we don't care about them anymore//

drop if notreat == 0

sum non_HDL_C - S_HDL_TG_pct if visit == 1 

//Put back in nice order//

xtset patient visit
order patient visit
sort patient visit

//Run log-linear model then exponentiate coefficients back to give change over time as a % change in geometric mean while adjusting for various covariates//

ds lognon_HDL_C - logS_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolites: `metabolite'"
 eststo: quietly mixed `metabolite' i.visit i.sex bl`metabolite' blbmi blage time|| patient:, reml 
 estadd expb 
   }
esttab _all using perc_change_over_time_addit_supplementary.csv, replace cells(expb(fmt(3)) se(fmt(3)) p(fmt(3))) gaps lines nostar keep(1.visit)
  eststo clear  

  

*******************************ASSOCIATION OF HBA1c WITH METABOLOMIC CHANGES ACROSS ADOLESCENCE IN UNTREATED PATIENTS*********************************

//Calculate residuals to use as exposures for increasing risk factors//

regress finalhba1c blhba1c
predict hba1c_res, res

regress finalbmi blbmi
predict bmi_res, res

zscore blhba1c finalhba1c blbmi finalbmi

//Final metabolite value based on baseline hba1c value//

ds z_non_HDL_C - z_S_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolites: `metabolite'"
  eststo: quietly mixed `metabolite' z_blhba1c blage i.sex blduration bl`metabolite' z_blbmi country || patient: if notreat == 1 
  }
esttab _all using hba1c_baseline_effect_supplementary.csv, replace cells(b(fmt(6)) se(fmt(6)) p(fmt(3))) gaps lines nostar keep(z_blhba1c)
  eststo clear
 
//Final metabolite value based on hba1c change across adolescence//

ds z_non_HDL_C - z_S_HDL_TG_pct
local metabolites `r(varlist)'
local i = 0
foreach metabolite of local metabolites {
  local i = `i' + 1
  di "Metabolites: `metabolite'"
  eststo: quietly mixed `metabolite' hba1c_res blage i.sex blduration bl`metabolite' bmi_res country || patient: if notreat == 1 
  }
esttab _all using hba1c_change_effect_supplementary.csv, replace cells(b(fmt(6)) se(fmt(6)) p(fmt(3))) gaps lines nostar keep(hba1c_res)
  eststo clear
  

