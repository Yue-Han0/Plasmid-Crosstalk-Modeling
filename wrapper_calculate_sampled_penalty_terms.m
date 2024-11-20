function [all_penalty_terms,penalty_labels,high_res_all_penalty_terms,high_res_all_penalty_terms_unscaled,high_res_penalty_term_labels,high_res_penalty_term_length] = wrapper_calculate_sampled_penalty_terms(sampled_params_selected,all_dosing_information,modelObj,experimental_grouped_data,kinetic_param_names,target_name_list)


    num_sample = size(sampled_params_selected,1); 
    num_init_cond = size(all_dosing_information,1); 

    % Evaluate model and calculate objective values 
    num_penalty_terms = 6; 
    all_penalty_terms = nan(num_sample * num_init_cond,num_penalty_terms); 
    high_res_num_penalty_terms = 32; 
    high_res_all_penalty_terms = nan(num_sample * num_init_cond,high_res_num_penalty_terms); 
    high_res_all_penalty_terms_unscaled = nan(num_sample * num_init_cond,high_res_num_penalty_terms); 
    % all_simulated_time = cell(num_sample * num_init_cond,1); 
    % all_simulated_data = cell(num_sample * num_init_cond,1); 
    problemObject = struct('Model',modelObj,'Data',experimental_grouped_data); 
    
    for sample_idx = 1:1%size(sampled_params_selected,1)
        nominal_sampled_params_single = sampled_params_selected(sample_idx,:); 
        % nominal_sampled_params_single = 10 .^ sampled_params_single; 
        
        % Sample initial conditions and create new dosing information 
        for init_sample_idx = 1:1%num_init_cond
            dosing_information = all_dosing_information{init_sample_idx}; 
    
            % Evaluate penalty values 
            [penalty_values,penalty_labels,~,~] = wrapper_calculate_obj_v2(modelObj,experimental_grouped_data,nominal_sampled_params_single,kinetic_param_names,dosing_information,target_name_list); 
            [high_res_scaled_penalty_vec,high_res_penalty_term_labels,high_res_penalty_term_length,high_res_unscaled_penalty_vec] = ...
                wrapper_calculate_obj_higher_resolution(nominal_sampled_params_single,problemObject,dosing_information,kinetic_param_names,target_name_list,'PE',7,4);
    
            all_penalty_terms((sample_idx - 1) * num_init_cond + init_sample_idx,:) = penalty_values; 
            high_res_all_penalty_terms((sample_idx - 1) * num_init_cond + init_sample_idx,:) = high_res_scaled_penalty_vec; 
            high_res_all_penalty_terms_unscaled((sample_idx - 1) * num_init_cond + init_sample_idx,:) = high_res_unscaled_penalty_vec; 
            % all_simulated_time{(sample_idx - 1) * num_init_cond + init_sample_idx,1} = simulated_time; 
            % all_simulated_data{(sample_idx - 1) * num_init_cond + init_sample_idx,1} = simulated_data; 
        end
    
    end

end