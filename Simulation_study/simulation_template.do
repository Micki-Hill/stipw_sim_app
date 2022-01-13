*****************************************
*			Simulation template			*
***************************************** 

* Note, macros have to be preceeded by \ 

* Date: 13/12/2021
* Version: v1.0


// Specify the details of the simulation study
local bootreps 500										// Number of bootstrap iterations
local adm 30											// Administrative censoring time
local delta 20											// Time at which difference in marginal RMST is calculcated										

mat cov = 	(1,0,0,0,0,0,0,0,0,0.05 \ ///
			 0,1,0,0,0,0,0,0,0,0.05 \ ///
			 0,0,1,0,0,0,0,0,0,0.05 \ ///
			 0,0,0,1,0,0,0,0,0,0.10 \ ///
			 0,0,0,0,1,0,0,0,0,0.20 \ ///
			 0,0,0,0,0,1,0,0,0,0.30 \ ///
			 0,0,0,0,0,0,1,0,0,0.10 \ ///
			 0,0,0,0,0,0,0,1,0,0.20 \ ///
			 0,0,0,0,0,0,0,0,1,0.30 \ ///
			 0.05, 0.05, 0.05, 0.10, 0.20, 0.30, 0.10, 0.20, 0.30, 1)
// I have put alpha 1 - 9 directly into the code
			 
			 
// Full list of each parameter that is varied		 
local samp_list 200 10000
local gamma_list 0.7 1.4
local lambda_list 0.15 0.003
local prev_list 0.1 0.25 0.5
local eff_list  0.5 2
local cens_list 0 20									// Mean of censoring distribution

// alpha_0 for each desired prevalence
local int_list -2.59 -1.33 0


// Bootstrap function

capture program drop usboot
program define usboot , rclass
	preserve
	bsample
	
	// Unstabilised weights
	qui logit a z1 z2 z3 z4 z5 z6 z7 z8 z9
	qui predict ps
	qui gen double uw = a/ps + (1-a)/(1-ps)
	qui stset y [pw = uw], failure(d)
	qui streg a, d(w)

	mat emat = e(b)
	local loghr = scalar(emat[1,colnumb(emat,"a")])
	local loglambda = scalar(emat[1,colnumb(emat,"_cons")])
	local gamma = scalar(exp(emat[1,colnumb(emat,"/ln_p")]))
	return scalar umhr = \`loghr'

	qui range ts 0 20 1001
	qui gen double usurv1 = exp(-exp(\`loglambda'+\`loghr')*ts^\`gamma')
	integ usurv1 ts
	local rmst1 = r(integral)
	return scalar urmst1 = \`rmst1'
	
	qui gen double usurv0 = exp(-exp(\`loglambda')*ts^\`gamma')
	integ usurv0 ts
	local rmst0 = r(integral)
	return scalar urmst0 = \`rmst0'
	
	return scalar udrmst = \`rmst1' - \`rmst0'
	

	// Stabilised weights
	qui sum a , meanonly
	qui local pr = r(mean)
	qui gen double sw = a*\`pr'/ps + (1-a)*(1-\`pr')/(1-ps)
	qui stset y [pw = sw], failure(d)
	qui streg a, d(w)

	mat emat = e(b)
	local loghr = scalar(emat[1,colnumb(emat,"a")])
	local loglambda = scalar(emat[1,colnumb(emat,"_cons")])
	local gamma = scalar(exp(emat[1,colnumb(emat,"/ln_p")]))
	return scalar smhr = \`loghr'

	qui gen double ssurv1 = exp(-exp(\`loglambda'+\`loghr')*ts^\`gamma')
	integ ssurv1 ts
	local rmst1 = r(integral)
	return scalar srmst1 = \`rmst1'
	
	qui gen double ssurv0 = exp(-exp(\`loglambda')*ts^\`gamma')
	integ ssurv0 ts
	local rmst0 = r(integral)
	return scalar srmst0 = \`rmst0'
	
	return scalar sdrmst = \`rmst1' - \`rmst0'
	
	restore
end		


// Function to repost variance so do not need to rerun stipw (saves time)
capture program drop robustvar
program define robustvar , eclass
	mat Vrob = e(V_robust)
	ereturn repost V = Vrob
end



// Simulation study

timer clear

// Post the random number state
capture postclose RNS
postfile RNS int(repno) int(gamma) int(prev) int(eff) int(cens)  int(samp) str2000(state1) str2000(state2) str2000(state3) using RandomNumberState_s\`scenario' , replace

// Post the estimations
capture postclose ests
postfile ests int(repno) int(gamma) int(prev) int(eff) int(cens) int(samp) str5(estimand) str1(weight) str6(method) ///
	float(aprev) float(cens_int) float(cens_adm) int(ntrt1) int(ntrt0) int(nevent1) int(nevent0) ///
	float(beta) float(se) int(error) int(converged) using estimates_s\`scenario' , replace

set seed \`seed' 
local n_samp : word count \`samp_list'


// Set the scenario
 
local gamma \`: word \`g' of \`gamma_list''
local lambda \`: word \`g' of \`lambda_list''
	
local prev \`: word \`p' of \`prev_list''
local int \`: word \`p' of \`int_list''
		
local eff \`: word \`e' of \`eff_list'' 
				
local cens \`: word \`c' of \`cens_list'' 


// Generate and analyse the data 

timer on 1
forvalues s = 1/\`n_samp' {
    local samp \`: word \`s' of \`samp_list''
	forvalues rep = 1/\`reps' {
		
		timer on 2
		// Get the random number state
		local state1 = substr(c(rngstate), 1,    2000)
		local state2 = substr(c(rngstate), 2001, 4000)
		local state3 = substr(c(rngstate), 4001, .)
		qui post RNS (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("\`state1'") ("\`state2'") ("\`state3'")
				
				
		// Data generation

		// Set the sample size
		clear
		qui set obs \`samp'

						
		// Create the covariates
		qui drawnorm z1 z2 z3 z4 z5 z6 z7 z8 z9 zu, cov(cov) double
							
							
		// Obtain variable U from Z_U
		qui gen double u = normal(zu)

						
		// Generate the treatment variable
		qui gen double p_a = 	(1+exp(-1*(\`int' + log(1.25)*z1 + log(1.50)*z2 + log(1.75)*z3 + ///
			log(1.25)*z4 + log(1.50)*z5 + log(1.75)*z6 + ///
			log(1.10)*z7 + log(1.10)*z8 + log(1.10)*z9)))^(-1)
		qui gen a = rbinomial(1,p_a)				
		
		
		// Monitor prevalences and number in each group
		qui sum a
		local aprev = r(mean)
		
		qui count if a == 1
		local ntrt1 = \`r(N)'			
		qui count if a == 0
		local ntrt0 = \`r(N)'					
				
				
		// Generate the outcome time variable (uses Weibull model)		
		qui gen double t = (-log(u)/(\`lambda' *exp(log(\`eff') *a) ))^(1/\`gamma')
				
				
		// Generate the event indicator variable
		qui gen double censt = rexponential(\`cens')
		qui gen double y = min(t, censt, \`adm')
		qui gen d = 0
		qui replace d = 1 if t < censt & t < \`adm'
			
		
		// Monitor censoring and number of events in each group
		qui count if d == 0 & y < \`adm'
		local cens_int = r(N)
		qui count if d == 0 & y == \`adm'
		local cens_adm = r(N)
		
		qui count if a == 1 & d == 1
		local nevent1 = \`r(N)'			
		qui count if a == 0 & d == 1
		local nevent0 = \`r(N)'
											
							
		// Analysis
						
		mat point = (.,.,.,.,.,.,.,.)							// For bootstrap point estimates
		foreach w in u s {
			
			qui stset y , failure(d)
			timer on 10
			cap qui stipw (logit a z1 z2 z3 z4 z5 z6 z7 z8 z9) , d(rp) df(1) noorthog ipw(\`w')
			timer off 10
							
			local rc_\`w' = _rc
			// Check for an error with stipw - M-estimation
			if (_rc >0) {
				foreach meth in "robust" "mest" {
					foreach estimand in "mhr" "rmst1" "rmst0" "drmst" {	
						qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("\`estimand'") ("\`w'") ("\`meth'") ///
							(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
							(.) (.) (1) (.)
					}
				}		
			}
			
			// If no error
			else {
				local converged = \`e(converged)'
				if ("\`w'" == "u") {							// For bootstrap point estimates
					local i = 0
				}
				else {
					local i = 4
				}
								
								
				// M-estimation
								
				// Marginal log hazard ratio (mhr)
				local rc_\`w'_mhr = 0
				qui mat emat = e(b)
				qui local \`w'mhr = scalar(emat[1,colnumb(emat,"a")])
				qui mat point[1,\`i'+1] = \`\`w'mhr'
				qui mat Vmest = e(V)
				qui local \`w'mhr_se = sqrt(scalar(Vmest[rownumb(Vmest,"a"), colnumb(Vmest,"a")]))
								
				// Marginal RMST in group 1 (rmst1)
				timer on 11
				cap qui predict \`w'rmst1 in 1, rmst at(a 1) stdp tmax(\`delta') 
				timer off 11
				local rc_\`w'_rmst1 = _rc
				if (\`rc_\`w'_rmst1' != 0) {
					qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("rmst1") ("\`w'") ("mest") ///
						(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
						(.) (.) (1) (.)
				}
				else {
				    qui local \`w'rmst1 = \`w'rmst1 in 1
					qui mat point[1,\`i'+2] = \`\`w'rmst1'
					qui local \`w'rmst1_se = \`w'rmst1_se in 1
				}
								
				// Marginal RMST in group 0 (rmst0)
				timer on 12
				cap qui predict \`w'rmst0 in 1, rmst at(a 0) stdp tmax(\`delta') 
				timer off 12
				local rc_\`w'_rmst0 = _rc
				if (\`rc_\`w'_rmst0' != 0) {
					qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("rmst0") ("\`w'") ("mest") ///
						(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
						(.) (.) (1) (.)
				}
				else {
					qui local \`w'rmst0 = \`w'rmst0 in 1
					qui mat point[1,\`i'+3] = \`\`w'rmst0'
					qui local \`w'rmst0_se = \`w'rmst0_se in 1	
				}
								
				// Marginal difference in RMST (drmst)
				timer on 13
				cap qui predictnl \`w'drmst = predict(rmst at(a 1) tmax(\`delta')) - predict(rmst at(a 0) tmax(\`delta')) in 1, se(\`w'drmst_se)
				timer off 13
				local rc_\`w'_drmst = _rc
				if (\`rc_\`w'_drmst' != 0) {
					qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("drmst") ("\`w'") ("mest") ///
						(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
						(.) (.) (1) (.)
				}
				else {
					qui local \`w'drmst = \`w'drmst in 1
					qui mat point[1,\`i'+4] = \`\`w'drmst'
					qui local \`w'drmst_se = \`w'drmst_se in 1
				}
				
				// Post the estimates from M-estimation for each of the estimands
				foreach estimand in mhr rmst0 rmst1 drmst {
				    if (\`rc_\`w'_\`estimand'' == 0) {
						qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("\`estimand'") ("\`w'") ("mest") ///
							(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
							(\`\`w'\`estimand'') (\`\`w'\`estimand'_se') (0) (\`converged')
					}
				}
				drop *rmst*
								
								
				// Robust
				robustvar										// Function defined at beginning to take the robust variance estimates from stipw, without having to re-run stipw
								
				// Marginal log hazard ratio (mhr)
				mat Vrob = e(V)
				qui local \`w'mhr_se = sqrt(scalar(Vrob[rownumb(Vrob,"a"), colnumb(Vrob,"a")]))
								
				// Marginal RMST in group 1 (rmst1)
				cap qui predict \`w'rmst1 in 1, rmst at(a 1) stdp tmax(\`delta') 
				local rc_\`w'_rmst1 = _rc
				if (\`rc_\`w'_rmst1' != 0) {
					qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("rmst1") ("\`w'") ("robust") ///
						(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
						(.) (.) (1) (.)
				}
				else {
					qui local \`w'rmst1_se = \`w'rmst1_se in 1
				}
								
				// Marginal RMST in group 0 (rmst0)
				cap qui predict \`w'rmst0 in 1, rmst at(a 0) stdp tmax(\`delta') 
				local rc_\`w'_rmst0 = _rc
				if (\`rc_\`w'_rmst0' != 0) {
					qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("rmst0") ("\`w'") ("robust") ///
						(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
						(.) (.) (1) (.)
				}
				else {
					qui local \`w'rmst0_se = \`w'rmst0_se in 1
				}
								
				// Marginal difference in RMST (drmst)
				cap qui predictnl \`w'drmst = predict(rmst at(a 1) tmax(\`delta')) - predict(rmst at(a 0) tmax(\`delta')) in 1, se(\`w'drmst_se)
				local rc_\`w'_drmst = _rc
				if (\`rc_\`w'_drmst' != 0) {
					qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("drmst") ("\`w'") ("robust") ///
						(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
						(.) (.) (1) (.)
				}
				else {
					qui local \`w'drmst_se = \`w'drmst_se in 1
				}
				
				// Post the estimates from robust for each of the estimands
				foreach estimand in mhr rmst0 rmst1 drmst {
				    if (\`rc_\`w'_\`estimand'' == 0) {
						qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("\`estimand'") ("\`w'") ("robust") ///
							(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
							(\`\`w'\`estimand'') (\`\`w'\`estimand'_se') (0) (\`converged')
					}
				}
				drop *rmst*
			}
		}
						
						
		// Bootstrapping	
		
		// Could not start bootstrap as at least one point estimate not available, error code 2				
		if (\`rc_u' != 0 | \`rc_s' != 0) {
		    foreach w in u s {
				foreach estimand in mhr rmst0 rmst1 drmst {
					qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("\`estimand'") ("\`w'") ("boot") ///
						(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
						(.) (.) (2) (.)	
				}
			}			
		}
		
		else {							
			timer on 20
			cap qui simulate umhr = r(umhr) urmst1 = r(urmst1) urmst0 = r(urmst0) udrmst = r(udrmst) ///
				smhr = r(smhr) srmst1 = r(srmst1) srmst0 = r(srmst0) sdrmst = r(sdrmst) , ///
				reps(\`bootreps')   : usboot
			timer off 20
			
			// Bootstrapping failed for at least one estimand for at least one weight, error code 1
			if (_rc > 0) {
			    foreach w in u s {
					foreach estimand in mhr rmst0 rmst1 drmst {
						qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("\`estimand'") ("\`w'") ("boot") ///
							(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
							(.) (.) (1) (.)	
					}
				}
			}
				
			else {
				timer on 21
				qui bstat, stat(point) n(\`samp')
				timer off 21
								
				mat Vboot = e(V)					
				foreach w in u s {
					// Post the estimates from bootstrapping for each of the estimands
					foreach estimand in mhr rmst0 rmst1 drmst {
						qui local \`w'\`estimand'_se = sqrt(scalar(Vboot[rownumb(Vboot,"\`w'\`estimand'"), colnumb(Vboot,"\`w'\`estimand'")]))
										
						qui post ests (\`rep') (\`g') (\`p') (\`e') (\`c') (\`samp') ("\`estimand'") ("\`w'") ("boot") ///
							(\`aprev') (\`cens_int') (\`cens_adm') (\`ntrt1') (\`ntrt0') (\`nevent1') (\`nevent0') ///
							(\`\`w'\`estimand'') (\`\`w'\`estimand'_se') (0) (\`converged')	
					}
				}	
			}
			
		} 
		timer off 2
	}
}
timer off 1


// Get last RNS
local rep = \`reps'+1
local state1 = substr(c(rngstate), 1,    2000)
local state2 = substr(c(rngstate), 2001, 4000)
local state3 = substr(c(rngstate), 4001, .)
qui post RNS (\`rep') (0) (0) (0) (0) (0) ("\`state1'") ("\`state2'") ("\`state3'")
postclose RNS

postclose ests

timer list