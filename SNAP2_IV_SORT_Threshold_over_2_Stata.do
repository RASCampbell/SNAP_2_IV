* Stata Do file to run the IV analysis for POMS estimate with postoperative CCU admission, for the SORT Threshold Sensitvity Analysis (>2%)

* This code runs a bivariate probit model to obtain the estimates as well as a 2SLS IV model to obtain indicators of model strength (Kleibergen-Paap Wald rank F statistic and Stock and Yogo critical values )

* The measures produced from the Bivariate probit model are: 
** Average Treatment Effect (AKA Average Partial Effect of "Treatment") = This is a marginal effect which describes the difference between the probability of our outcome occuring, given ccu admisison, and the probability of our outcome occuring, given ward admission. Thus, ATE is the difference between the marginal probability of outcome success given treatment "success" and the marginal probability of outcome success given treatment "failure". This can be descibed in lay terms as "on average, patients admitted to icu are X% less/more likely to have POMS/mortality". This is done by holding the other covariates at a constant point and changing the treatment variable.

** Probit coefficient & CI = The coefficient, on the other hand, provides the estimated effect of a specific explanatory variable on the probability of an event occurring in one equation, without considering causality. It represents the association between the explanatory variable and the outcome, but it does not necessarily imply a causal relationship. A positive coefficient means that an increase in the predictor leads to an increase in the predicted probability. A negative coefficient means that an increase in the predictor leads to a decrease in the predicted probability. However, interpretation of the coefficients in probit regression is not as straightforward as the interpretations of coefficients in linear regression or logit regression.  The increase in probability attributed to a one-unit increase in a given predictor is dependent both on the values of the other predictors and the starting value of the given predictors.

** Odds ratio and CI = Generally the magnitudes of coefficent estimates should not be compared accross different kinds of models as they use different scale factors. However, as a rough rule of thumb, multiplying probit coefficients by 1.6 makes them roughtly comparable to logit estimates (i.e. log odds), which is what we have done here to aid interpretation. 

clear all


* Import data using dta files
use patients_labelled6_imputed_sort2_over, clear

* Summarise the data * Check that all variables are correctly coded
describe

summarize  CCUCapacityTimeofSurgery  i.POMS i.POMS_factor i.icu_adm_factor   i.S01Gender   Age_scaled  i.age_p10 i.weekend_factor i.S02OperativeUrgency   i.S02PlannedProcSeverity   i.asa   i.PreopLOS_factor  i.combined_DHx   i.combined_PMHx   i.diabetes_factor sysBP_scaled  HR_scaled i.Anaesthetist_grade_combind_fctr   i.Surgeon_grade_combined_factor i.abnormal_bloodtst_mpt_nrml_fctr i.S04BloodLoss 

* create icu_adm dummy variable 
generate icu_adm_dummy = (icu_adm_factor==2)

* Create list of observed confounding variables to adjust the model for
global xlist ///
		i.S01Gender ///
		i.age_p10 ///
		i.weekend_factor ///
		i.S02OperativeUrgency   ///
		i.S02PlannedProcSeverity ///
		i.asa ///
		i.PreopLOS_factor ///
		i.combined_DHx ///
		i.combined_PMHx ///
		i.diabetes_factor ///
		sysBP_scaled ///
		HR_scaled   ///
		i.Anaesthetist_grade_combind_fctr   ///
		i.Surgeon_grade_combined_factor  ///
		i.abnormal_bloodtst_mpt_nrml_fctr  ///
		i.S04BloodLoss

* Run bivariate probit model 
biprobit (outcome: POMS = i.icu_adm_dummy $xlist) (rxselect: icu_adm_dummy = $xlist CCUCapacityTimeofSurgery), nolog  vce(cluster SiteCode)

est store iv_biprobit

* ATE (aka marginal effect)
margins, dydx(icu_adm_dummy) predict(pmarg1) force
* ATE are manually extracted from here and inputed into R Script for plotting

* save bivariate probit results in an interpretable table
tempfile estimates_file
parmest, ///
	label list(parm label estimate min* max* p) ///
	stars(0.05 0.01 0.001) ///
	format(estimate min* max* %9.2f p %9.3f) ///
	saving(`estimates_file', replace)

use `estimates_file', clear
save $table_name.dta, replace

*  ============
*  = Latexify =
*  ============

use $table_name.dta, clear
d
list eq parm est min95 max95 stderr
gen eq_order = 0 if eq == "rxselect"



replace eq_order = 1 if eq == "outcome"
gen eq_label  = "Selection model" if eq == "rxselect"
replace eq_label  = "Outcome model" if eq == "outcome"
drop if parm == "_cons"

* convert to odds ratios via factor 1.6
gen or_coef = 1.6 * est
gen or_stderr   = 1.6 * stderr
gen or_est    = exp(or_coef)
gen or_min95    = exp(or_coef - invnormal(0.975) * or_stderr)
gen or_max95    = exp(or_coef + invnormal(0.975) * or_stderr)

list eq parm est min95 max95 or_est or_min95 or_max95

* Relabel variables (especialy ones with multiple levels)
cap drop parm_label
gen parm_label = ""
replace parm_label = "CCU adm False" if parm == "0b.icu_adm_dummy"
replace parm_label = "CCU adm True" if parm == "1.icu_adm_dummy"
replace parm_label = "Female" if parm == "1b.S01Gender"
replace parm_label = "Male" if parm == "2.S01Gender"
replace parm_label = "18 to 30" if parm == "1b.age_p10"
replace parm_label = "31 to 40" if parm == "2.age_p10"
replace parm_label = "41 to 50" if parm == "3.age_p10"
replace parm_label = "51 to 60" if parm == "4.age_p10"
replace parm_label = "61 to 70" if parm == "5.age_p10"
replace parm_label = "71 to 80" if parm == "6.age_p10"
replace parm_label = "80+" if parm == "7.age_p10"
replace parm_label = "Weekday" if parm == "1b.weekend_factor"
replace parm_label = "Weekend" if parm == "2.weekend_factor"
replace parm_label = "Procedure count = 1" if parm == "1b.procedure_count_factor"
replace parm_label = "Procedure count = 2" if parm == "2.procedure_count_factor"
replace parm_label = "Procedure count > 2" if parm == "3.procedure_count_factor"
replace parm_label = "Elective" if parm == "1b.S02OperativeUrgency"
replace parm_label = "Expedited" if parm == "2.S02OperativeUrgency"
replace parm_label = "Urgent" if parm == "3.S02OperativeUrgency"
replace parm_label = "Immediate" if parm == "4.S02OperativeUrgency"
replace parm_label = "Minimal" if parm == "1b.S02PlannedProcSeverity"
replace parm_label = "Intermediate" if parm == "2.S02PlannedProcSeverity"
replace parm_label = "Major" if parm == "3.S02PlannedProcSeverity"
replace parm_label = "Very Major" if parm == "4.S02PlannedProcSeverity"
replace parm_label = "Complex" if parm == "5.S02PlannedProcSeverity"
replace parm_label = "ASA I" if parm == "1b.asa"
replace parm_label = "ASA II" if parm == "2.asa"
replace parm_label = "ASA III" if parm == "3.asa"
replace parm_label = "ASA IV/V" if parm == "4.asa"
replace parm_label = "GA False" if parm == "1b.S04AnaestheticTechniqueGeneral"
replace parm_label = "GA True" if parm == "2.S04AnaestheticTechniqueGeneral"
replace parm_label = "Pre-op Level of Support = 0" if parm == "1b.S02LevelOfSupport"
replace parm_label = "Pre-op Level of Support = 1" if parm == "2.S02LevelOfSupport"
replace parm_label = "Pre-op Length of Stay = 0" if parm == "1b.PreopLOS_factor"
replace parm_label = "Pre-op Length of Stay = 1" if parm == "2.PreopLOS_factor"
replace parm_label = "Pre-op Length of Stay = 2 to 6" if parm == "3.PreopLOS_factor"
replace parm_label = "Pre-op Length of Stay = 7 to 20" if parm == "4.PreopLOS_factor"
replace parm_label = "Pre-op Length of Stay = 21 to 251" if parm == "5.PreopLOS_factor"
replace parm_label = "CXR cardiomegaly" if parm == "1b.combined_radiology_findings"
replace parm_label = "CXR consolidation" if parm == "2.combined_radiology_findings"
replace parm_label = "CXR not done" if parm == "3.combined_radiology_findings"
replace parm_label = "CXR normal" if parm == "4.combined_radiology_findings"
replace parm_label = "CXR other abnormality" if parm == "5.combined_radiology_findings"
replace parm_label = "Drug history False" if parm == "1b.combined_DHx"
replace parm_label = "Drug history True" if parm == "2.combined_DHx"
replace parm_label = "PMH False" if parm == "1b.combined_PMHx"
replace parm_label = "PMH True" if parm == "2.combined_PMHx"
replace parm_label = "Not diabetic" if parm == "1b.diabetes_factor"
replace parm_label = "Type 2 diabetes - diet controlled" if parm == "2.diabetes_factor"
replace parm_label = "Type 2 diabetes - oral medication" if parm == "3.diabetes_factor"
replace parm_label = "Type 2 diabetes - insulin" if parm == "4.diabetes_factor"
replace parm_label = "Type 1 diabetes" if parm == "5.diabetes_factor"
replace parm_label = "GCS 8 or less" if parm == "1b.GCS_category_factor"
replace parm_label = "GCS 9-11" if parm == "2.GCS_category_factor"
replace parm_label = "GCS 12-14" if parm == "3.GCS_category_factor"
replace parm_label = "GCS 15" if parm == "4.GCS_category_factor"
replace parm_label = "Pre-op BP scaled" if parm == "sysBP_scaled"
replace parm_label = "Pre-op HR scaled" if parm == "HR_scaled"
replace parm_label = "Pre-op dyspnoea level = None" if parm == "1b.dyspnoea_factor"
replace parm_label = "Pre-op dyspnoea level = On exertion" if parm == "2.dyspnoea_factor"
replace parm_label = "Pre-op dyspnoea level = Limiting activity" if parm == "3.dyspnoea_factor"
replace parm_label = "Pre-op dyspnoea level = At rest" if parm == "4.dyspnoea_factor"
replace parm_label = "Night surgery = False" if parm == "1b.night_surgery_factor"
replace parm_label = "Night surgery = True" if parm == "2.night_surgery_factor"
replace parm_label = "Most senior anaesthetist grade = Consultant" if parm == "1b.Anaesthetist_grade_combind_fctr"
replace parm_label = "Most senior anaesthetist grade = Non-consultant" if parm == "2.Anaesthetist_grade_combind_fctr"
replace parm_label = "Most senior surgeon grade = Consultant" if parm == "1b.Surgeon_grade_combined_factor"
replace parm_label = "Most senior surgeon grade = Non-consultant" if parm == "2.Surgeon_grade_combined_factor"
replace parm_label = "Pre-assessment = No" if parm == "1b.pre_assessment_include_na"
replace parm_label = "Pre-assessment = Yes" if parm == "2.pre_assessment_include_na"
replace parm_label = "Pre-assessment = Not applicable" if parm == "3.pre_assessment_include_na"
replace parm_label = "Abnormal blood tests = False" if parm == "1b.abnormal_bloodtst_mpt_nrml_fctr"
replace parm_label = "Abnormal blood tests = True" if parm == "2.abnormal_bloodtst_mpt_nrml_fctr"
replace parm_label = "Critical unexpected events perioperatively = False" if parm == "1b.S04CriticalUnxpctdEvntsPrprtvly"
replace parm_label = "Critical unexpected events perioperatively = True" if parm == "2.S04CriticalUnxpctdEvntsPrprtvly"
replace parm_label = "Intra-op blood loss = 0-100ml" if parm == "1b.S04BloodLoss"
replace parm_label = "Intra-op blood loss = 101-500ml" if parm == "2.S04BloodLoss"
replace parm_label = "Intra-op blood loss = 501-999ml" if parm == "3.S04BloodLoss"
replace parm_label = "Intra-op blood loss = >1000ml" if parm == "4.S04BloodLoss"
replace parm_label = "Available CCU beds at time of surgery" if parm == "CCUCapacityTimeofSurgery"

list eq parm parm_label est min95 max95 p stars or_est or_min95 or_max95

* export excel sheet with results
export excel using "CCU_IV_sort_2+_pars_mv.xlsx", sheet("raw") sheetreplace firstrow(varlabels)


* Run 2SLS IV analysis
clear all

* Import data using dta files
use patients_labelled6_imputed_sort2_over, clear

* Summarise the data * Check that all variables are correctly coded
describe

summarize  CCUCapacityTimeofSurgery  i.POMS i.POMS_factor i.icu_adm_factor   i.S01Gender   Age_scaled  i.age_p10 i.weekend_factor i.S02OperativeUrgency   i.S02PlannedProcSeverity   i.asa   i.PreopLOS_factor  i.combined_DHx   i.combined_PMHx   i.diabetes_factor sysBP_scaled  HR_scaled i.Anaesthetist_grade_combind_fctr   i.Surgeon_grade_combined_factor i.abnormal_bloodtst_mpt_nrml_fctr i.S04BloodLoss 

* create icu_adm dummy variable 
generate icu_adm_dummy = (icu_adm_factor==2)

* Create list of observed confounding variables to adjust the model for
global xlist ///
		i.S01Gender ///
		i.age_p10 ///
		i.weekend_factor ///
		i.S02OperativeUrgency   ///
		i.S02PlannedProcSeverity ///
		i.asa ///
		i.PreopLOS_factor ///
		i.combined_DHx ///
		i.combined_PMHx ///
		i.diabetes_factor ///
		sysBP_scaled ///
		HR_scaled   ///
		i.Anaesthetist_grade_combind_fctr   ///
		i.Surgeon_grade_combined_factor  ///
		i.abnormal_bloodtst_mpt_nrml_fctr  ///
		i.S04BloodLoss


** 2SLS 
ivreg2 POMS	 ///
		$xlist ///
		(i.icu_adm_dummy = CCUCapacityTimeofSurgery) ///
		, ///
		cluster(SiteCode) first rf endog(i.icu_adm_dummy)

* Save results into Excel sheet
tempfile estimates_file
parmest, ///
	label list(parm label estimate min* max* p) ///
	escal(widstat jp) ///
	stars(0.05 0.01 0.001) ///
	format(estimate min* max* %9.4f p %9.4f) ///
	saving(`estimates_file', replace)

use `estimates_file', clear
save $table_name.dta, replace

use $table_name.dta, clear

sdecode estimate, format(%9.3fc) gen(est)
sdecode min95, format(%9.3fc) replace
sdecode max95, format(%9.3fc) replace
sdecode p, format(%9.3fc) replace
replace p = "<0.01" if p == "0.00"
replace p = "<0.001" if p == "0.000"
gen est_ci95 = "(" + min95 + " -- " + max95 + ")" if !missing(min95, max95)

* export excel sheet with results
export excel using "2SLS_CCU_IV_sort_2+_pars_mv.xlsx", sheet("raw") sheetreplace firstrow(varlabels)


