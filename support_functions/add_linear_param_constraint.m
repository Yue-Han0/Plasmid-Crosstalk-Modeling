function [inequality_A,inequality_b] = add_linear_param_constraint(param_oi_list,base_inequality_matrix,base_inequality_rhs,Estimated,existing_inequality_A,existing_inequality_B)
    % Convert heuristic/known linear relationship/constraints for estimated
    % parameters into the constraint form in optimization functions 

    % INPUT
        % param_oi_list: a cell array with names of parameters invovled in the constraint 
        % base_inequality_matrix: a # constraints * length(param_oi_list) matrix representing A_eq 
        % base_inequality_rhs: a # constraints * 1 vector containing RHS of the linear constraints 
        % Estimated: an array of estimatedInfo object in Simbiology containing all parameters to be estimated 

    % OUTPUT 
        % inequality_A: a # constriants * # estimated parameters array for A input in optimization solver such as fmincon 
        % inequality_b: a # constriants * # estimated parameters array for b input in optimization solver such as fmincon 

    num_estimated = length(Estimated); 
    estimated_params_name_list = {Estimated.Name}; 
    
        % Find location in estimated parameters 
    param_idx_list = nan(length(param_oi_list),1); 
    for param_oi_idx = 1:length(param_oi_list)
        param_oi_name = param_oi_list{param_oi_idx};
        param_idx_list(param_oi_idx) = find(strcmp(param_oi_name,estimated_params_name_list)); 
    end
    
        % Convert base unequality matrix
    temp_inequality_A = zeros(size(base_inequality_matrix,1),num_estimated); 
    for inequality_idx = 1:length(base_inequality_rhs)
        temp_inequality_A(inequality_idx,param_idx_list) = base_inequality_matrix(inequality_idx,:); 
    end

    if ~exist('existing_inequality_A','var')
        inequality_A = temp_inequality_A; 
    else
        inequality_A = [existing_inequality_A;temp_inequality_A]; 
    end
    if ~exist('existing_inequality_B','var')
        inequality_b = base_inequality_rhs; 
    else
        inequality_b = [existing_inequality_B;base_inequality_rhs];
    end

end