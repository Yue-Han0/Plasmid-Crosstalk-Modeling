function PACE_wrapper_toxin_RNAP
    currentpath = pwd;
    addpath(genpath(sprintf('%s/plasmid_crosstalk_config_files',currentpath))); 
    addpath(genpath(sprintf('%s/txtlsim_buildacell',currentpath))); 
    addpath(genpath(sprintf('%s',currentpath)));
    rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath));

    group_description = {'PE_T7_strong_no_empty'};
    param_info_path = 'benchmark_model/parameters_test_PE.xlsx';
    
    weight = []; 
    
    problemObject_RNAP = setProblemObject_PE_var1_toxRNAP(group_description,param_info_path,weight); 
    
    wrapper_fit_to_data_w_probObject(problemObject_RNAP,'RNAP_toxin_mech')

end