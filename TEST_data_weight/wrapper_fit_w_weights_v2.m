function wrapper_fit_w_weights_v2

    group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
    'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
    'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
    'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
    param_info_path = 'parameters_v2.xlsx';
    problemObject_weights = setProblemObject_v3_normalized(group_description,[],param_info_path);
    save_path_id_weight = 'w_weight'; 
    override_result_path = 'param_est_run_save/20240829_param_est_run1207_w_weight';


    wrapper_fit_to_data_w_probObject_v2(problemObject_weights,save_path_id_weight,override_result_path)


end