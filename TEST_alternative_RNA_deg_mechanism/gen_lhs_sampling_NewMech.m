function gen_lhs_sampling_NewMech(mechanism)

    
    % Sample using the updated parameter ranges 
    group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
        'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
        'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
        'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
    param_info_path = 'param_info/parameters_v2_expanded_newMech.xlsx';
    if strcmp(mechanism,'TestBindingSites')
        problemObject = setProblemObject_v3_TestBindingSite(group_description,[],param_info_path); 
    elseif strcmp(mechanism,'TestPrototype')
        problemObject = setProblemObject_v3_TestPrototype(group_description,[],param_info_path); 
    end
    dosing_information = create_dosing_info_from_problemObject(problemObject); 
    modelObj = problemObject.Model; 
    experimental_grouped_data = problemObject.Data; 

    num_kinetic_params = length(problemObject.Estimated); 
    
    kinetic_param_names = {problemObject.Estimated.Name};
    kinetic_param_bounds = [problemObject.Estimated.Bounds]; 
    kinetic_param_bounds = reshape(kinetic_param_bounds,[2,num_kinetic_params]); 
    kinetic_param_bounds = kinetic_param_bounds'; 
    
    log_param_bounds = log10(kinetic_param_bounds); 
    
    % Generate sampled kinetic parameter values 
    num_sample = 1000000; % 1e+06
    lhsSamples = lhsdesign(num_sample,num_kinetic_params); 
    sampled_params = nan(num_sample,num_kinetic_params); 
    for i = 1:num_kinetic_params
        lowerBound = log_param_bounds(i, 1);
        upperBound = log_param_bounds(i, 2);
        sampled_params(:, i) = lowerBound + (upperBound - lowerBound) * lhsSamples(:, i);
    end
    
    nominal_sampled_params = 10 .^ sampled_params;
    
        % Define parameters to be constrained/filtered 
    kinetic_param_name_oi_list = {'TXTL_PT7_RNAPbound_R','TXTL_PT773_RNAPbound_R',...
        'TXTL_PJ23119_RNAPbound_R','TXTL_PJ23105_RNAPbound_R'}; 
        % Find their indices in params list 
    kinetic_param_name_oi_idx_list = nan(length(kinetic_param_name_oi_list),1); 
    for kinetic_param_oi_idx = 1:length(kinetic_param_name_oi_idx_list)
        kinetic_param_name_oi = kinetic_param_name_oi_list{kinetic_param_oi_idx}; 
        kinetic_param_name_oi_idx_list(kinetic_param_oi_idx) = find(strcmp(kinetic_param_name_oi,kinetic_param_names)); 
    end
    
        % Encode a matrix to represent relationships among parameters
            % Each row represents a criteria and each column represent
            % parameters involved 
    criteria_matrix = [1 -1  0  0;
                       0  1  0  -1;
                       0  0  1  -1;]; % Assumes < 0  
    satisfied_idx_list = []; % Initialize 
    for criterion_idx = 1:size(criteria_matrix,1)
        criteria_vec = criteria_matrix(criterion_idx,:); 
        smaller_param_idxa = find(criteria_vec == 1); 
        larger_param_idxa = find(criteria_vec == -1); 
        row_satisfied_idx_list = nominal_sampled_params(:,kinetic_param_name_oi_idx_list(smaller_param_idxa)) < nominal_sampled_params(:,kinetic_param_name_oi_idx_list(larger_param_idxa));
        if isempty(satisfied_idx_list)
            satisfied_idx_list = row_satisfied_idx_list; 
        else
            satisfied_idx_list = satisfied_idx_list & row_satisfied_idx_list; 
        end
    
    end
    final_satisfied_idx_list = find(satisfied_idx_list);
    
    satisfied_sampled_params = nominal_sampled_params(final_satisfied_idx_list,:); 
    
    save(sprintf('test_save_files/20241007_lhs_sampled_params_%s.mat',mechanism),'satisfied_sampled_params','problemObject','dosing_information'); 

end