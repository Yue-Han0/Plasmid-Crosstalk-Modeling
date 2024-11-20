function param_summary = get_all_estimated_params(all_fitResults)
    % Get a # parameters * # runs matrix of estimated parameters from
    % fitting results 

    % Preallocation 
    num_run = length(all_fitResults); 
    num_param = height(all_fitResults{1,1}.ParameterEstimates); 
    param_summary = nan(num_param,num_run);

    % Assign fitted parameters 
    for run_idx = 1:num_run
        fitResult = all_fitResults{run_idx,1}; 
        param_summary(:,run_idx) = [fitResult.ParameterEstimates.Estimate];
    end



end