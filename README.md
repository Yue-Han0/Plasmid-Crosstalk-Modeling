# Plasmid-Crosstalk-Modeling
 
## Overview
Code repository for *manuscript A Mathematical Model for Cell-Free Transcription-Translation for Plasmid Crosstalk*    
Authors: Yue Han, Alexandra T. Patterson, Fernanda Piorino, Mark P. Styczynski    
    
## 1. Data Processing     
All data and code related are stored in the *Data* folder, which contains the following directories:     
**raw_data**: includes raw data from plate reader     
**data_preprocessing**: converts raw data into matlab-friendly format. Results saved in processed data folder    
**processed_data**: described above      
**calibration_curve**: standard curves for 3WJdB and sfGFP    
**data_structures**: organized data tables containing processed experimental data    
**plasmid_crosstalk_config_files**: Excel files defining biology parts to be used in txtlsim toolbox    

## 2. Model Generation     
Construct Simbiology model object containing experimental data, mechanistic model, and parameter information    
**param_info**: contains information on the initial guesses, lower & upper bound of model parameters         
**setProblemObject_v3.m**: construct the described model object        

## 3. Model Calibration 
**wrapper_fit_qualitative_characteristics_w_probObject.m**: model calibration to minimize sum of surrogate terms     
**wrapper_fit_to_data_w_probObject_v2.m**: model calibration to minimize deviation      

## 4. Surrogate Terms for Alternative Model Evaluation 
**surrogate_terms**: functions to calculate surrogate terms based on multiple sets of time-course data    
**wrapper_calculate_obj_higher_resolution.m**: output all surrogate terms based on model and parameters     
**wrapper_calculate_obj_higher_resolution_from_data.m**: output all surrogate terms based on time-course data     

## 5. Parameter Sampling 
**run_param_and_init_cond_sampling_v2.m**: calculate surrogate terms for lhs sampled parameters 

## 6. Other experiments 
### 6.1 Alternative RNA degradation mechanism     
### 6.2 Model calibration with normalized data     
### 6.3 Alternative toxin mechanism

## 7. Results & Visualization 
**result_files**: contains surrogate term values with estimated parameters or sampled parameters 
**plot_functions**：contains functions for visualization
