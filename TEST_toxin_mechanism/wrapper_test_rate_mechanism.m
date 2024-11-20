function wrapper_test_rate_mechanism


    group_description = {'PE_T7_strong_no_empty'};
    
    param_info_path = 'benchmark_model/parameters_test_PE.xlsx';

    save_path_id_list = {'TXgen_Ribo','TLgen_Ribo','TLgen_RNAP'}; 
    
    
    weight = []; 
    
    for test_idx = 1:length(save_path_id_list)
        save_path_id = save_path_id_list{test_idx};
        switch test_idx
            case 1
                problemObject = setProblemObject_PE_var2_TXgen_Ribo(group_description,param_info_path,weight); 
            case 2 
                problemObject = setProblemObject_PE_var3_TLgen_Ribo(group_description,param_info_path,weight); 
            case 3
                problemObject = setProblemObject_PE_var4_TLgen_RNAP(group_description,param_info_path,weight); 
        end
        wrapper_fit_to_data_w_probObject(problemObject,save_path_id);
    end

end