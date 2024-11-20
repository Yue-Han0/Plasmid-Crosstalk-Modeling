
function [dev,penalty] = calculate_baseline_data_dev(experimental_time,experimental_data,simulated_time,simulated_data)
    
    % Extract baseline dataset index list 
    baseline_idx_list = [];
    for baseline_idx = 1:28:length(simulated_data)-1
        baseline_idx_list = [baseline_idx_list baseline_idx:baseline_idx+6];
    end
    baseline_experimental_time = experimental_time(baseline_idx_list);
    baseline_experimental_data = experimental_data(baseline_idx_list);
    baseline_simulated_time = simulated_time(baseline_idx_list);
    baseline_simulated_data = simulated_data(baseline_idx_list);

    % Extrapolate simulated data to timeVec of the experimental data
    dev = 0; 
    for dataset_idx = 1:length(baseline_experimental_data)
        baseline_experimental_time_single = baseline_experimental_time{dataset_idx}; 
        baseline_experimental_data_single = baseline_experimental_data{dataset_idx}; 
        baseline_simulated_time_single = baseline_simulated_time{dataset_idx}; 
        baseline_simulated_data_single = baseline_simulated_data{dataset_idx}; 

        extrap_baseline_simulated_data_single = interp1(baseline_simulated_time_single(1:10:end),baseline_simulated_data_single(1:10:end,1),...
            baseline_experimental_time_single,'linear','extrap'); 
        dev = dev + sum((baseline_experimental_data_single - extrap_baseline_simulated_data_single(:,1)).^2); 
    end

    % calculate penalty to regularize decreasing expression at higher
    % reporter plasmid concentrations 
        % Going to semi-hardcode this 
    T7_strong_15 = baseline_simulated_data{6}; 
    T7_strong_30 = baseline_simulated_data{7}; 
    T7_weak_15 = baseline_simulated_data{13}; 
    T7_weak_30 =  baseline_simulated_data{14}; 
    sigma70_strong_10 = baseline_simulated_data{19};
    sigma70_strong_30 = baseline_simulated_data{21}; 

    penalty = max(T7_strong_30(:,1)) - max(T7_strong_15(:,1)) + max(T7_weak_30(:,1)) - max(T7_weak_15(:,1))...
        + max(sigma70_strong_30(:,1)) - max(sigma70_strong_10(:,1)); 

end