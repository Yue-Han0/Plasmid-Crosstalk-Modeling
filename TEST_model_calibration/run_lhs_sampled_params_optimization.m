function run_lhs_sampled_params_optimization(run_idx)
    
    % Add path 
    currentpath = pwd;
    addpath(genpath(sprintf('%s',currentpath)));
    rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath));

    % Initialize sampled parameters 
    init_param_file = load('test_save_files/20240822_lhs_sampled_params_optimization.mat'); 
    satisfied_sampled_params = init_param_file.satisfied_sampled_params;
    problemObject = init_param_file.problemObject;

    % Select a manageable pool of parameters for optimization 
    num_param_set_per_run = 100; 
    sampled_params_selected = satisfied_sampled_params((run_idx - 1) * num_param_set_per_run + 1:run_idx * num_param_set_per_run,:);

    % Set save_path_id and run parameters estimation 
    save_path_id = sprintf('param_est_run_save/20240822_lhs_sampled_params_opt_qual_obj_run%d.mat',run_idx); 
    wrapper_fit_to_data_w_probObject_params(problemObject,sampled_params_selected,save_path_id)


end