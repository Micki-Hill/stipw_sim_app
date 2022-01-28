README file for stipw paper
# Variance estimation in inverse probability weighted parametric survival models

In the manuscript titled _Variance estimation in inverse probability weighted parametric survival models_ by Hill et al (in draft), a framework for a closed-form variance estimator for inverse probability (IP) weighted parametric survival models is proposed. The method can be implemented using a newly developed Stata command _stipw_, which can be obtained from SSC or from the [`stipw`] (link) repository. This repository contains the Stata code for the corresponding simulation study, which evaluated the performance of the variance estimator, and an illustration of its use on data from an AIDS clinical trial.


## Simulation study

The methodology was evaluated using a simulation study. The three aims of the simulation study were:
* To evaluate the performance of the M-estimation variance estimator (for both stabilised and unstabilised weights)
* To compare whether stabilised or unstabilised weights result in a better M-estimation variance estimator
* To compare the computational time of the M-estimation and bootstrap variance estimators

The simulation study was performed for 24 scenarios. The simulation study focused on large sample properties, where a sample size of 10000 was used. As an exploratory analysis, a sample size of 200 was also investigated. Each scenario (including both sample sizes) was ran as a separate batch on the High Performance Computing cluster at the University of Leicester. 

[`create_simulation_files.do`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/create_simulation_files.do) shows how the simulation _.do_ files were created. This was performed twice: once for the iteration sample size calculation (1000 repetitions) and once for the main simulation study (7700 repetitions). The file uses the [`simulation_template.do`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/simulation_template.do) to produce the 24 simulation _.do_ files (one for each scenario). An example is shown with [`simulation_s1.do`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/simulation_s1.do) for scenario 1 for the main simulation study.

The 24 datasets (one for each scenario) from the main simulation study are merged together in [`analyse_simulation_study.do`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/analyse_simulation_study.do) and the difficult entries discussed in Section 6 in the manuscript are removed. The combined dataset is then transformed into wide format and split into four .csv files (split by sample size and gamma value). This was done so that the files were small enough to be uploaded onto GitHub. [`analyse_simulation_study.do`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/analyse_simulation_study.do) includes code to rebuild the combined dataset. To do this, the user needs to download, unzip and save the following .csv files:
* [`estimates_g1_s200.7z`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/estimates_g1_s200.7z)
* [`estimates_g2_s200.7z`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/estimates_g2_s200.7z)
* [`estimates_g1_s10000.7z`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/estimates_g1_s10000.7z)
* [`estimates_g2_s10000.7z`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/estimates_g2_s10000.7z)

[`analyse_simulation_study.do`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Simulation_study/analyse_simulation_study.do) performs the simulation study analysis in Section 4.2 of the manuscript. This includes:
* Merging the simulation datasets
* Removing the difficult entries discussed in Section 6
* Exporting the dataset into .csv files
* Rebuilding the combined dataset
* Calculating the performance measures with MCSE
* Checking the MCSE is within the desired tolerance as discussed in Supplementary Material S1
* Creating an example table of the results as given in Supplementary Material S4
* Creating an example graph of the results as given in Section 4.2 (and Supplementary Material S3)


## Application

The methodology was illustrated with the AIDS Clinical Trial Group Study 175 (ACTG175) dataset. The randomised clinical trial compared monotherapy with zidovudine or didanosine with combination therapy with zidovudine and didanosine or zidovudine and zalcatbine in adults with HIV type I. For more information, please see Hammer et al. 

The ACTG175 dataset can be obtained from the [speff2trial](https://CRAN.R-project.org/package=speff2trial)  package in R and consists of 2139 individuals and 17 possible confounders. First, the dataset needs to be exported as a .csv file. This can be done in R with the following code:

```
install.packages(“speff2trial”)
library(speff2trial)
write.csv(ACTG175, file=”actg175.csv”, row.names=F)
```

[`actg175_application.do`](https://github.com/Micki-Hill/stipw_sim_app/blob/main/Application/actg175_application.do) performs the analysis in Chapter 5 of the manuscript. This includes:
* Initial formatting of the data
* Fitting the propensity score model and creating and investigating the weights
* Checking the standardised differences of covariates before and after weighting
* Summarising the survival data
* Model selection for the IP weighted outcome survival model
* Point estimation with and without weighting
* Variance estimation for stabilised and unstabilised weights: robust, M-estimation and bootstrapping (M=500 and M=10000)

#### References
* Hammer SM, Katzenstein DA, Hughes MD, et al. A trial comparing nucleoside monotherapy with combination therapy in HIV-infected adults with CD4 cell counts from 200 to 500 per cubic millimeter. N Engl J Med. 1996;335(15):1081-1090.
