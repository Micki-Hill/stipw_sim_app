* Aim: To analyse the ACTG175 dataset as described in Chapter 5 in the manuscript

* Date: 13/12/2021
* Version: v1.0

* Steps performed:
* 	1. Initial formatting of the data
* 	2. Fitting the propensity score model and creating and investigating the weights
* 	3. Checking the standardised differences of covariates before and after weighting
* 	4. Summarising the survival data
* 	5. Model selection for the IP weighted outcome survival model
* 	6. Point estimation with and without weighting
* 	7. Variance estimation for stabilised and unstabilised weights: 
*			(i)   Robust
*			(ii)  M-estimation
*			(iii) Bootstrapping (M=500 and M=10000)

* Trt: monotherapy (arms 0 and 3) vs combination therapy (arms 1 and 2)
* Outcome: years until the first occurrence of: 
*		(i)   a decline in CD4 T cell count of at least 50 
*		(ii)  an event indicating progression to AIDS 
*		(iii) death
* Confounders: all 15 baseline confounders except zprior (as was singular) and 
*			   str2 (antihist) as collinearity with strat


*****************************************
*			Data formatting				*
*****************************************

// First export the ACTG175 dataset from the speff2trial package in R, see the README file

import delimited using "actg175.csv", clear


rename v14 antihist																
// Does not allow str2 as variable name

label define arms 0 "Zidovudine" 1 "Zidovudine and didanosine" 2 "Zidovudine and zalcitabine" 3 "Didanosine"
label values arms arms

label define race 0 "White" 1 "Non-White"
label values race race

label define gender 0 "Female" 1 "Male"
label values gender gender

// Change time-to-event variable to be in years
gen double years = days / 365.25

// Recode treatment to binary to compare monotherapy and combination therapy
gen trt = 0 if arms == 0 | arms == 3
replace trt = 1 if arms == 1 | arms == 2
label define trt 0 "Monotherapy" 1 "Combined"
label values trt trt

// Create dummy variables for strat
tab strat, gen(strat)

save actg175, replace


*****************************************
* 		Propensity score & weights		*
*****************************************

// Number treated and untreated
tab trt

// Propensity score model
logit trt race gender age wtkg homo drugs hemo symptom z30 preanti strat2 strat3 ///
	cd40 cd80 oprior karnof
predict ps	

// Unstabilised weights
gen double uw = trt/ps + (1-trt)/(1-ps)
sum uw, d

// Stabilised weights
count if trt == 1
local nt = `r(N)'
count 
local n = `r(N)'
gen double sw = `nt'/`n' * trt/ps  + (1-`nt'/`n')*(1-trt)/(1-ps)
sum sw, d


*****************************************
*		Standardised differences		*
*****************************************

// Before PS (split into continuous and factor variables)
foreach cvar in age wtkg karnof preanti cd40 cd80 {
	pstest `cvar', raw t(trt)
}
foreach fvar in race gender homo drugs hemo symptom oprior z30 strat {
	pstest i.`fvar', raw t(trt)
}
// Look at B statistic: Max of 6.4 (0.064)


* After PS: unstabilised weights
foreach cvar in age wtkg karnof preanti cd40 cd80 {
	pstest `cvar', t(trt) mw(uw)
}
foreach fvar in race gender homo drugs hemo symptom oprior z30 strat {
	pstest i.`fvar', t(trt) mw(uw)
}
// Look at B statistic: Max of 3.4 (0.034)


* After PS: stabilised weights
foreach cvar in age wtkg karnof preanti cd40 cd80 {
	pstest `cvar', t(trt) mw(sw)
}
foreach fvar in race gender homo drugs hemo symptom oprior z30 strat {
	pstest i.`fvar', t(trt) mw(sw)
}
// Look at B statistic: Max of 3.4 (0.034)


*****************************************
*		Summary survival data	 		*
*****************************************

// Person years
stset years, failure(cens)

// Maximum follow-up
sum years

// How many events in each group
tab trt cens , row


*****************************************
*			Model selection		 		*
*****************************************

// Unstabilised weights

stset years [pw = uw], failure(cens)

// Standard distributions
foreach dist in exp weib gomp logl logn {
	streg trt, d(`dist')
	estimates store u_`dist'
}

// RP models with 1-5 degrees of freedom
forvalues i = 1/5 {
	// trt as time-independent
	stpm2 trt, scale(hazard) df(`i')
	estimates store u_rp`i'
	
	// trt as time-dependent with 1-5 degrees of freedom
	forvalues j = 1/5 {
		stpm2 trt, scale(hazard) df(`i') tvc(trt) dftvc(`j') 
		estimates store u_rp`i'_`j'		
	}
}

count
estimates stat u*, n(`r(N)')


// Stabilised weights

stset years [pw = sw], failure(cens)

// Standard distributions
foreach dist in exp weib gomp logl logn {
	streg trt, d(`dist')
	estimates store s_`dist'
}

// RP models with 1-5 degrees of freedom
forvalues i = 1/5 {
	// trt as time-independent
	stpm2 trt, scale(hazard) df(`i')
	estimates store s_rp`i'
	
	// trt as time-dependent with 1-5 degrees of freedom
	forvalues j = 1/5 {
		stpm2 trt, scale(hazard) df(`i') tvc(trt) dftvc(`j') 
		estimates store s_rp`i'_`j'		
	}
}

count
estimates stat s*, n(`r(N)')

// A Royston-Parmar (stpm2) model with 2 degrees of freedom and 1 degree of freedom for trt was chosen
// 5th lowest AIC and 1st lowest BIC for unstabilised, 2nd lowest AIC and 4th lowest BIC with stabilised


*****************************************
*			Point estimation			*
*****************************************

// Calculate difference in marginal restricted mean survival time (RMST) at time 3

// Without weighting (marginal as a randomised clinical trial)
stset years , failure(cens)
stpm2 trt, scale(hazard) df(2) tvc(trt) dftvc(1) eform
predictnl drmst = predict(rmst at(trt 1) tmax(3)) - predict(rmst at(trt 0) tmax(3)) in 1, se(drmst_se)

// Unstabilised weights
stset years [pw = uw], failure(cens)
stpm2 trt, scale(hazard) df(2) tvc(trt) dftvc(1) eform
predictnl udrmst = predict(rmst at(trt 1) tmax(3)) - predict(rmst at(trt 0) tmax(3)) in 1

// Stabilised weights
stset years [pw = sw], failure(cens)
stpm2 trt, scale(hazard) df(2) tvc(trt) dftvc(1) eform
predictnl sdrmst = predict(rmst at(trt 1) tmax(3)) - predict(rmst at(trt 0) tmax(3)) in 1

list *drmst* in 1
drop *drmst*


*****************************************
*	Variance estimation: Unstabilised	*
*****************************************

// Install stipw, if not already done
*ssc install stipw

timer clear
stset years, failure(cens)


// Robust
timer on 1
stipw (logit trt race gender age wtkg homo drugs hemo symptom z30 preanti strat2 strat3 ///
	cd40 cd80 oprior karnof) , ///
	dist(rp) df(2) dftvc(1) vce(robust) ipw(u) noheader
predictnl udrmst_rob = predict(rmst at(trt 1) tmax(3)) - predict(rmst at(trt 0) tmax(3)) in 1, se(udrmst_rob_se)
timer off 1
list udrmst_rob_se in 1


// M-estimation
timer on 2
stipw (logit trt age wtkg hemo homo drugs karnof oprior z30 preanti race gender ///
	strat2 strat3 symptom cd40 cd80) , ///
	dist(rp) df(2) dftvc(1) vce(mestimation) ipw(u) noheader
predictnl udrmst_mest = predict(rmst at(trt 1) tmax(3)) - predict(rmst at(trt 0) tmax(3)) in 1, se(udrmst_mest_se)
timer off 2
list udrmst_mest_se in 1


// Bootstrapping

// Make point estimate for bootstrap
local drmst = udrmst_rob in 1
mat point = (`drmst')

// Bootstrap program
capture program drop uboot
program define uboot , rclass
	use actg175, replace
	bsample
	
	// Weights
	qui logit trt race gender age wtkg homo drugs hemo symptom z30 preanti strat2 strat3 ///
		cd40 cd80 oprior karnof
	qui predict ps
	qui gen double uw = trt/ps + (1-trt)/(1-ps)
	
	// Outcome model
	qui stset years [pw = uw], failure(cens)
	qui stpm2 trt, scale(hazard) df(2) tvc(trt) dftvc(1)

	// RMST 1
	predict rmst1 in 1, rmst at(trt 1) tmax(3)
	local rmst1 = rmst1 in 1
	
	// RMST 0
	predict rmst0 in 1, rmst at(trt 0) tmax(3)
	local rmst0 = rmst0 in 1
	
	// Difference in RMST
	return scalar udrmst = `rmst1' - `rmst0'
end		

// Bootstrap: 500 reps
timer on 3
set seed 37493639
simulate udrmst = r(udrmst) , reps(500)  : uboot	
bstat, stat(point) n(2139)
timer off 3

// Bootstrap: 10000 reps
timer on 4
set seed 45027302	
simulate udrmst = r(udrmst) , reps(10000)  : uboot 	
bstat, stat(point) n(2139)
timer off 4


// Check timings - may be different to text as depends on computer speed
timer list


*****************************************
*	Variance estimation: Stabilised		*
*****************************************

timer clear
use actg175, replace
stset years, failure(cens)


// Robust
timer on 1
stipw (logit trt race gender age wtkg homo drugs hemo symptom z30 preanti strat2 strat3 ///
	cd40 cd80 oprior karnof) , ///
	dist(rp) df(2) dftvc(1) vce(robust) ipw(s) noheader
predictnl sdrmst_rob = predict(rmst at(trt 1) tmax(3)) - predict(rmst at(trt 0) tmax(3)) in 1, se(sdrmst_rob_se)
timer off 1
list sdrmst_rob_se in 1


// M-estimation
timer on 2
stipw (logit trt age wtkg hemo homo drugs karnof oprior z30 preanti race gender ///
	strat2 strat3 symptom cd40 cd80) , ///
	dist(rp) df(2) dftvc(1) vce(mestimation) ipw(s) noheader
predictnl sdrmst_mest = predict(rmst at(trt 1) tmax(3)) - predict(rmst at(trt 0) tmax(3)) in 1, se(sdrmst_mest_se)
timer off 2
list sdrmst_mest_se in 1

	
// Bootstrapping

// Make point estimate for bootstrap
local drmst = sdrmst_rob in 1
mat point = (`drmst')

// Bootstrap program
capture program drop sboot
program define sboot , rclass
	use actg175, replace
	bsample
	
	* Weights
	qui logit trt race gender age wtkg homo drugs hemo symptom z30 preanti strat2 strat3 ///
		cd40 cd80 oprior karnof
	qui predict ps
	qui sum trt , meanonly
	qui local pr = r(mean)
	qui gen double sw = trt*`pr'/ps + (1-trt)*(1-`pr')/(1-ps)
	
	* Outcome model
	qui stset years [pw = sw], failure(cens)
	qui stpm2 trt, scale(hazard) df(2) tvc(trt) dftvc(1)

	* RMST 1
	predict rmst1 in 1, rmst at(trt 1) tmax(3)
	local rmst1 = rmst1 in 1
	
	* RMST 0
	predict rmst0 in 1, rmst at(trt 0) tmax(3)
	local rmst0 = rmst0 in 1
	
	* Difference in RMST
	return scalar sdrmst = `rmst1' - `rmst0'
end		

// Bootstrap: 500 reps
timer on 3
set seed 836154383
simulate  sdrmst = r(sdrmst) , reps(500)  : sboot	
bstat, stat(point) n(2139)
timer off 3

// Bootstrap: 10000 reps
timer on 4
set seed 37856298	
simulate sdrmst = r(sdrmst) , reps(10000)  : sboot 	
bstat, stat(point) n(2139)
timer off 4


// Check timings - may be different to text as depends on computer speed
timer list