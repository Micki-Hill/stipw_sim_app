* Aim: To analyse the results of the main simulation study as described in Section 4.2 in the manuscript

* Date: 28/01/2022
* Version: v1.1

* Steps performed:
* 	1. Merge the simulation datasets
* 	2. Remove the difficult entries discussed in Section 6
*	3. Export the dataset into .csv files
*	4. Rebuild the estimates dataset
* 	5. Calculate the performance measures with MCSE
* 	6. Check the MCSE is within the desired tolerance as discussed in Supplementary Material S1
* 	7. Create an example table of the results as given in Supplementary Material S4
* 	8. Create an example graph of the results as given in Section 4.2 (and Supplementary Material S3)

* Note that steps 1-3 are for reference. Steps 4-8 can be ran by the user, once they have
* downloaded, unzipped and saved the 4 csv files from GitHub.


/*

*************************************
*		Merge the datasets			*
*************************************

// This section appends the 24 datasets, one created for each scenario.

use estimates_s1, replace
forvalues i = 2/24 {
	append using estimates_s`i'
}

qui gen methodn = 1 if weight == "u" & method == "robust" 
qui replace methodn = 2 if weight == "u" & method == "boot" 
qui replace methodn = 3 if weight == "u" & method == "mest" 
qui replace methodn = 4 if weight == "s" & method == "robust" 
qui replace methodn = 5 if weight == "s" & method == "boot"
qui replace methodn = 6 if weight == "s" & method == "mest"
qui label define methodn 1 "U: Robust" 2 "U: Boot" 3 "U: M-est" 4 "S: Robust" 5 "S: Boot" 6 "S: M-est"
qui label values methodn methodn


*************************************
*		Removal of difficulties		*
*************************************

// There were no convergence issues for the large sample size.
// This section removes the difficult results for the small sample size as discussed in Section 6 in the manuscript. 


// Manual input of results for 2 datasets (small sample size)

// There were 2 incidents where unstabilised weights gave an error in the simulation study (for all methods), but produced results outside of the 
// simulation study (the dataset was reproduced using the seed). The simulation study results were updated with the correct results (code not shown).

// This affected rep 226 and 1656 where gamma = 1.4, treatment prevalence = 10%, hazard ratio = 0.5 and there was intermittent censoring. 


// Remove difficulties 1: Error caused by few individuals in the treatment group leading to extreme bootstrap samples that did not converge
drop if error == 1
// Removed 168 observations (21 datasets x 4 estimands x 2 bootstrap methods).

// Remove difficulties 2: No events in at least one treatment group
drop if nevent1 == 0 | nevent0 == 0
// Should remove 17688 (737 datasets x 4 estimands x 6 methods).
// Removed 17672 as 16 (2 datasets x 4 estimands x 2 bootstrap) were already removed as they had error = 1.

// Remove difficulties 3: Considerably large model based standard errors caused by few events in one treatment group leading to extreme bootstrap samples
drop if se > 30
// Removed 6 observations.


save estimates, replace


*************************************
*		Export as csv files			*
*************************************

// This section transforms the data into wide format and exports it as four .csv files (one for each gamma and sample size).
// This is done so that the files are small enough to be uploaded onto GitHub.

use estimates, replace

drop method weight																// methodn combines these
reshape wide beta se error converged , ///
	i(repno gamma prev eff cens samp estimand) j(methodn) 

order repno gamma prev eff cens samp estimand aprev cens_int cens_adm ntrt1 ntrt0 nevent1 nevent0

export delimited "estimates_g1_s200.csv" if samp == 200 & gamma == 1, replace
export delimited "estimates_g2_s200.csv" if samp == 200 & gamma == 2, replace
export delimited "estimates_g1_s10000.csv" if samp == 10000 & gamma == 1, replace
export delimited "estimates_g2_s10000.csv" if samp == 10000 & gamma == 2, replace

*/


*************************************
*		Rebuild estimates dataset	*
*************************************

// This section imports and rebuilds the estimates dataset for analysis.
// Please download, unzip and save each of the 4 csv datasets from Github.

// Import the csv files
import delimited "estimates_g1_s200.csv", clear
save estimates_g1_s200, replace

import delimited "estimates_g2_s200.csv", clear
save estimates_g2_s200, replace

import delimited "estimates_g1_s10000.csv", clear
save estimates_g1_s10000, replace

import delimited "estimates_g2_s10000.csv", clear
save estimates_g2_s10000, replace

// Append the csv files
use estimates_g1_s200, replace
append using estimates_g2_s200 
append using estimates_g1_s10000 
append using estimates_g2_s10000

// Reshape the csv files
reshape long beta se error converged , ///
	i(repno gamma prev eff cens samp estimand) j(methodn)
qui label define methodn 1 "U: Robust" 2 "U: Boot" 3 "U: M-est" 4 "S: Robust" 5 "S: Boot" 6 "S: M-est"
qui label values methodn methodn

// Drop the estimates that had previously been removed (now they are just missing)
drop if beta == .

save estimates, replace

		
*************************************
*   Calculate performance measures	*
*************************************

// This section calculates the relative error of ModSE with MCSE and CIs.

// Calculate EmpSE, ModSE^2 (Model based variance) and the SD of ModSE
use estimates, replace
gen double var = se^2
collapse (sd) empse = beta (mean) modvar = var (sd) sdmodvar = var (count) nrep = repno, by(estimand gamma prev eff cens samp methodn)

// Calculate ModSE and the relative error
gen double modse = sqrt(modvar)
gen double relerror = 100*(modse/empse-1)

// Calculate the MCSE
gen double varmodvar = sdmodvar^2
gen double mcse = 100*modse/empse * sqrt(varmodvar/(4*nrep*modse^4) + 1/(2*(nrep-1)))
// nrep is used (rather than 7700) as some observations were removed for the small sample and so nrep was less than 7700 in these cases.

// Calculate the confidence intervals for the relative error
gen double lci = relerror - 1.96*mcse
gen double uci = relerror + 1.96*mcse

save estimates_summary, replace


*************************************
*		Check the MCSE				*
*************************************

// This section checks the maximum MCSE across the estimand, method, weight and scenario combinations is < 1 percentage point.
// This corresponds to Supplementary Material S1.

use estimates_summary, replace
bys estimand: sum mcse if samp == 10000
// They are within the threshold of 1 percentage point.


*************************************
*			Results: Tables			*
*************************************

// This section formats the data into the style of the tables given in Supplementary Material S4.

use estimates_summary, replace

keep empse relerror mcse estimand gamma prev eff cens methodn samp

rename relerror po_value1
rename mcse po_value2
reshape long po_value, i(estimand gamma prev eff cens samp methodn empse) j(po)

reshape wide empse po_value , i(estimand gamma prev eff cens samp po) j(methodn)
label define po 1 "Point" 2 "MCSE"
label values po po

forvalues i = 1/6 {
	tostring empse`i', replace force format(%9.3f)
	tostring po_value`i', replace force format(%9.2f)
	replace po_value`i' = "(" + po_value`i' + ")" if po == 2
}
order empse1 po_value1 po_value2 po_value3 empse4 po_value4 po_value5 po_value6
sort samp gamma prev eff cens samp po
	
// The following gives the table for the marginal log hazard ratio, sample size of 10000
br gamma prev eff cens empse1 - po_value6 if estimand == "mhr" & samp == 10000


*************************************
*			Results: Graphs			*
*************************************

// This section provides example code to create a results graph, this code is for Figure 2.
// The results graphs are discussed in Section 4.2 in the manuscript.

set scheme s1mono
use estimates_summary, replace
qui gen ref = 0

gen scenario = 1 if eff == 1 & cens == 1
replace scenario = 2 if eff == 1 & cens == 2
replace scenario = 3 if eff == 2 & cens == 1
replace scenario = 4 if eff == 2 & cens == 2

replace scenario = scenario - 0.15 if methodn == 1 | methodn == 4
replace scenario = scenario + 0.15 if methodn == 3 | methodn == 6

qui label define prev 1 "Prev: 10%" 2 "Prev: 25%" 3 "Prev: 50%"
qui label values prev prev

qui label define gamma 1 "Gamma: 0.7" 2 "Gamma: 1.4"
qui label values gamma gamma

qui label define scenario 1 "0.5, NIC" 2 "0.5, IC" 3 "2, NIC" 4 "2, IC" 
qui label values scenario scenario

local samp 10000
local estimand mhr					

twoway (rspike uci lci scenario if estimand == "`estimand'" & samp == `samp' & methodn == 4, by(gamma prev, ///
	title("Marginal log hazard ratio, Stabilised weights, N = `samp'") note("")) lcolor(maroon)) ///
	(rspike uci lci scenario if estimand == "`estimand'" & samp == `samp' & methodn == 5, by(gamma prev) lcolor(navy)) ///
	(rspike uci lci scenario if estimand == "`estimand'" & samp == `samp' & methodn == 6, by(gamma prev) lcolor(dkgreen)) ///
	(scatter relerror scenario if estimand == "`estimand'" & samp == `samp' & methodn == 4, by(gamma prev) msymbol(smcircle)  msize(small) mcolor(maroon)) ///
	(scatter relerror scenario if estimand == "`estimand'" & samp == `samp' & methodn == 5, by(gamma prev) msymbol(smcircle)  msize(small) mcolor(navy)) ///
	(scatter relerror scenario if estimand == "`estimand'" & samp == `samp' & methodn == 6, by(gamma prev) msymbol(smcircle)  msize(small) mcolor(dkgreen)) ///
	, ///
	xtitle("Scenario", size(small))  ytitle("Relative % Error", size(small))  ///
	xlabel(1(1)4, valuelabel labsize(small)) ylabel(,labsize(small)) graphregion(margin(vsmall)) ///
	legend(order(4 "Robust" 5 "Bootstrap" 6 "M-estimation") rows(1) size(small)) yline(0, lcolor(black))
