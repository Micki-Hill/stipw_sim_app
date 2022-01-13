* Aim: To create the 24 simulation files for the main simulation study

* Date: 13/12/2021
* Version: v1.0

// Ensure stipw is installed
ssc install stipw

clear all

// Specify how many repetitions
*local reps 1000					// Iteration sample size calculation
local reps 7700						// Main simulation study

/*
// Seeds for 1000 run
local seed_list ///
	68974653 ///
	92831467 ///
	69967318 ///
	58266073 ///
	9822963 ///
	65929101 ///
	70364889 ///
	56090126 ///
	86491290 ///
	92904620 ///
	99863536 ///
	46156778 ///
	64649368 ///
	99484752 ///
	34727798 ///
	60356018 ///
	34348778 ///
	35748674 ///
	78631317 ///
	14177532 ///
	28509040 ///
	13538605 ///
	56200428 ///
	15548638
*/

// Seeds for 7700 run
local seed_list ///
	8285889 ///
	442154 ///
	8630377 ///
	4526046 ///
	7720400 ///
	5861199 ///
	4227766 ///
	2729307 ///
	8053644 ///
	3060019 ///
	1190997 ///
	7247310 ///
	6964866 ///
	9119345 ///
	6795634 ///
	3549416 ///
	7389700 ///
	1874017 ///
	3146128 ///
	1375693 ///
	6537739 ///
	2701319 ///
	8998394 ///
	5734232
	
local n_seed : word count `seed_list'

// Create the simulation files for each of the 24 scenarios
local i = 1
forvalues g = 1/2 {
	forvalues p = 1/3 {
		forvalues e = 1/2 {
			forvalues c = 1/2 {
			
				file open sim using simulation_s`i'.do, write text replace
				local seed `: word `i' of `seed_list''
				
				// Write the details of the scenario to the file
				file write sim "local reps `reps'" _newline
				file write sim "local seed `seed'" _newline
				file write sim "local scenario `i'" _newline
				file write sim "local g `g'" _newline							// Gamma and lambda DGM parameters
				file write sim "local p `p'" _newline							// Treatment prevalence
				file write sim "local e `e'" _newline							// Treatment effect
				file write sim "local c `c'" _newline _newline					// Censoring

				// Copy in the simulation template
				file open template using simulation_template.do, read		
				file read template line
				while r(eof)==0 {
					file write sim `"`line'"' _newline
					file read template line
				}
				file close template
				
				file close sim
				local i = `i' + 1
			}
		}
	}	
}