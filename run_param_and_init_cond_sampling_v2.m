function run_param_and_init_cond_sampling_v2(run_idx)
    
    currentpath = pwd;
    addpath(genpath(sprintf('%s',currentpath)));
    rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath));
    
    init_param_file = load('test_save_files/20240819_lhs_sampled_params.mat'); 
    satisfied_sampled_params = init_param_file.satisfied_sampled_params;
    problemObject = init_param_file.problemObject_updated;
    dosing_information = init_param_file.dosing_information; 

    num_param_set_per_run = 100; 

    sampled_params_selected = satisfied_sampled_params((run_idx - 1) * num_param_set_per_run + 1:run_idx * num_param_set_per_run,:);

    high_res_all_penalty_terms = nan(num_param_set_per_run,32);
    high_res_all_penalty_terms_unscaled = nan(num_param_set_per_run,32);
    mode = 'PE';
    num_conc = 7;
    num_promotor = 4; 

    for i = 1:num_param_set_per_run 
        sampled_params_single = sampled_params_selected(i,:);
        [scaled_penalty_vec,penalty_term_labels,penalty_term_length,unscaled_penalty_vec] = wrapper_calculate_obj_higher_resolution(sampled_params_single,...
            problemObject,dosing_information,[],[],mode,num_conc,num_promotor);
        high_res_all_penalty_terms(i,:) = scaled_penalty_vec; 
        high_res_all_penalty_terms_unscaled(i,:) = unscaled_penalty_vec; 
    end

    
    save(sprintf('param_est_run_save/20241029_param_sampling_constrained_run%d.mat',run_idx),'high_res_all_penalty_terms','high_res_all_penalty_terms_unscaled',...
    'sampled_params_selected','dosing_information');

end