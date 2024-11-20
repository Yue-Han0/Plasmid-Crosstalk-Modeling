function lhs_sample(problemObject,kinetic_param_bounds,num_sample,save_path)
    % Sample parameters in the LHS space and save 

    % INPUT
        % problemObject: a fitProblem struct in Simbiology
        % kinetic_param_bounds: a # params * 2 matrix representing LB and
        % UB of the parameters 
        % num_sample: a double for number of samples in LHS 
        % save_path: a string for path/to/save
    
    %% Sanity Check
    num_kinetic_params = length(problemObject.Estimated); 
    if ~isequal(num_kinetic_params,size(kinetic_param_bounds,1))
        error('Inconsistent number of parameters')
    end
    
    %% Generate sampled kinetic parameter values 
    log_param_bounds = log10(kinetic_param_bounds); 
    lhsSamples = lhsdesign(num_sample,num_kinetic_params); 
    sampled_params = nan(num_sample,num_kinetic_params); 
    for i = 1:num_kinetic_params
        lowerBound = log_param_bounds(i, 1);
        upperBound = log_param_bounds(i, 2);
        sampled_params(:, i) = lowerBound + (upperBound - lowerBound) * lhsSamples(:, i);
    end
    
    nominal_sampled_params = 10 .^ sampled_params;
    
    %% Define parameters to be constrained/filtered 
    kinetic_param_name_oi_list = {'TXTL_RNAdeg_kc_sfGFP','TXTL_RNAdeg_kc_kanR','TXTL_RNAdeg_kc',...
        'TXTL_RNAdeg_R_sfGFP','TXTL_RNAdeg_R','TXTL_PT7_RNAPbound_R','TXTL_PT773_RNAPbound_R',...
        'TXTL_PJ23119_RNAPbound_R','TXTL_PJ23105_RNAPbound_R'}; 
    kinetic_param_names = {problemObject.Estimated.Name};
        % Find their indices in params list 
    kinetic_param_name_oi_idx_list = nan(length(kinetic_param_name_oi_list),1); 
    for kinetic_param_oi_idx = 1:length(kinetic_param_name_oi_idx_list)
        kinetic_param_name_oi = kinetic_param_name_oi_list{kinetic_param_oi_idx}; 
        kinetic_param_name_oi_idx_list(kinetic_param_oi_idx) = find(strcmp(kinetic_param_name_oi,kinetic_param_names)); 
    end
    
        % Encode a matrix to represent relationships among parameters
            % Each row represents a criteria and each column represent
            % parameters involved 
    criteria_matrix = [1  0  -1  0  0  0  0  0  0 ; 
                       0  1  -1  0  0  0  0  0  0 ;
                       0  0   0 -1  1  0  0  0  0 ;
                       0  0   0  0  0  1 -1  0  0 ;
                       0  0   0  0  0  0  1  0  -1;
                       0  0   0  0  0  0  0  1  -1;]; % Assumes < 0  
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
    
    save(save_path,'satisfied_sampled_params','problemObject'); 




end