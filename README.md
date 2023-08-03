Unveiling the Unobservable Causal Inference on Multiple Derived Outcomes
================



## Data

### Description

We use the fMRI data publicly available at the Autism Brain Imaging Data Exchange I.
79 ASD patients from the New York University Langone Medical Center are used in the Case Study. 
For each patient, there are 175 scans of brain overtime. A brain scan at one time point 
contains the brain signals of 116 regions-of-interests (ROI) extracted according to the 
automated anatomical labeling atlas. Covariates information including age, sex, handedness 
and FIQ are also used in the data analysis. The code used to analyze the data is provided.

### Process to request the data

The data are publicly available by request. 

The process to request the data is shown in readme file under the `data` folder.


## Code

### Version of primary software used

R version 4.1.3

### R packages needed

openxlsx_4.2.4

MASS_7.3-55

iterators_1.0.14

doParallel_1.0.17

foreach_1.5.2  

glmnet_4.1-4

doRNG_1.8.2  


### Code Files 

1.Functions.R is the code of the functions for simulation, including the function for data
generating, function for proposed method, function for BH procedure.

2.Simulation.R is the code for simulation of the proposed method on correlation. 

3.Simulation_regression.R is the code for simulation of the proposed method on regression parameters 
of lasso and debiased-lasso. 

4.plot.R uses simulation results in Simulation.R and Simulation_regression.R to produce Figure 1, 2 and S1.

5.RealCaseStudy.R is the code for application. This code produces table 1 and the brain connection in 
Figure S2.


## Reproducibility workflow

#### To reproduce the Case Study result:

First load in the COV_NYU.xlsx and every NYU_XXXXXXX_rois_aal.1D documents
(this step is also inlcuded in RealCaseStudy.R) 

Then use the code in RealCaseStudy.R to reproduce the analysis results. 

In RealCaseStudy.R, `Connection_diff1' provides the Nerwork in Table1 and Figure S2, 
`Tau_sig' provides the `Estimated Effect' in Table1, `CI_low' and `CI_high' provides 
the 95% CI in Table1, `Ave-trt' and `Ace-cl' provide Ave-trt and Ace-cl in Table1.

#### To reproduce the Simulation results: 

First load in the functions in Functions.R. 

Then use the code in Simulation.R  to reproduce the simulation results for the subject level correlations. In Simulation.R,
`Simulation_block_proposed' is the result of proposed method under block diagonal setting, 
`Simulation_block_BH' is the result of BH method under block diagonal setting, 
`Simulation_off_proposed' is the result of proposed method under off diagonal setting, 
`Simulation_off_BH' is the result of BH method under off diagonal setting. 
The simulation results are further organized in simulation-result.xlsx to produce Figure 1 and Figure 2.

Next, use the code in simulation_regression.R to reproduce the simulation results for the regression coeffcients. In Simulation_regression.R, `result' provides the analysis result used to construct Figure S1. The simulation results are further organized in simulation-result.xlsx to produce Figure S1.

Finally, use plot.R and the simulation results in simulation-result.xlsx to produce the figures. 
In plot.R, `plot_block' produces Figure 1, `plot_off' produces Figure 2, and `plot_lasso' produces Figure S1.



