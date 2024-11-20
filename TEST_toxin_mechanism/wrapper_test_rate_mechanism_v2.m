function wrapper_test_rate_mechanism_v2


    group_description = {'PE_T7_strong_no_empty'};
    
    param_info_path = 'benchmark_model/parameters_v2.xlsx';

    save_path_id_list = 'TXgen_Ribo_v2'; 
    
    problemObject = setProblemObject_PE_var2_TXgen_Ribo_v2(group_description,group_number,param_info_path); 
    wrapper_fit_to_data_w_probObject(problemObject,save_path_id);

end