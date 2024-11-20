function [scaled_penalty_vec,penalty_term_labels,penalty_term_length,unscaled_penalty_vec] = wrapper_calculate_obj_higher_resolution_from_data(simulated_time,simulated_data,...
    problemObject,mode,num_conc,num_promotor)

    % Extract information from problemObject
    experimental_grouped_data = problemObject.Data;
    species_name_list = {problemObject.Model.Species.Name};

        % Calculate crosstalk ratios from simulated data 
    simulated_crosstalk_ratios = calculate_crosstalk_ratio_v2(simulated_time,simulated_data,num_conc,num_promotor,mode); 

    [experimental_time,experimental_data] = process_crosstalk_data_from_source(experimental_grouped_data,'grouped_data'); 
    experimental_crosstalk_ratios = calculate_crosstalk_ratio_v2(experimental_time,experimental_data,num_conc,num_promotor,mode); 

    %% Calculate terms in penalty_vec 
        % penalty_vec = (1) baseline_data_dev (2) crosstalk_ratio_dev (3)
        % T7strong toxin (4) T7weak toxin (5) sigma70strong toxin (6)
        % sigma70weak no toxin (7) positive crosstalk (8) negative
        % crosstalk (9) promotor strength (apply to positive crosstalk
        % only) (10) other penalties

        %% (1) Baseline Data Deviation 
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
    baseline_dev = 0; 
    for dataset_idx = 1:length(baseline_experimental_data)
        baseline_experimental_time_single = baseline_experimental_time{dataset_idx}; 
        baseline_experimental_data_single = baseline_experimental_data{dataset_idx}; 
        baseline_simulated_time_single = baseline_simulated_time{dataset_idx}; 
        baseline_simulated_data_single = baseline_simulated_data{dataset_idx}; 

        extrap_baseline_simulated_data_single = interp1(baseline_simulated_time_single(1:10:end),baseline_simulated_data_single(1:10:end,1),...
            baseline_experimental_time_single,'linear','extrap'); 
        baseline_dev = baseline_dev + sum((baseline_experimental_data_single - extrap_baseline_simulated_data_single(:,1)).^2); 
    end

        %% (2) Crosstalk Ratio Deviation 
        % Calculate the deviation of crosstalk ratios 
    crosstalk_ratio_dev = 0; 
    for crosstalk_ratio_row_idx = 1:size(experimental_crosstalk_ratios,1)
        for crosstalk_ratio_col_idx = 1:size(experimental_crosstalk_ratios,2)
            experimental_crosstalk_ratio = experimental_crosstalk_ratios{crosstalk_ratio_row_idx,crosstalk_ratio_col_idx};
            simulated_crosstalk_ratio = simulated_crosstalk_ratios{crosstalk_ratio_row_idx,crosstalk_ratio_col_idx}; 
            crosstalk_ratio_dev = crosstalk_ratio_dev + sum((simulated_crosstalk_ratio - experimental_crosstalk_ratio).^2); 
        end
    end

        %% (3) T7 Strong Baseline Penalty 

    T7_strong_15 = baseline_simulated_data{6}; 
    T7_strong_30 = baseline_simulated_data{7}; 
    T7_strong_toxin = max(T7_strong_30(:,1)) - max(T7_strong_15(:,1)); 

        %% (4) T7 Weak Baseline Penalty
    T7_weak_15 = baseline_simulated_data{13}; 
    T7_weak_30 =  baseline_simulated_data{14}; 
    T7_weak_toxin = max(T7_weak_30(:,1)) - max(T7_weak_15(:,1)); 

        %% (5) sigma70 strong baseline penalty
    sigma70_strong_10 = baseline_simulated_data{19};
    sigma70_strong_30 = baseline_simulated_data{21}; 
    sigma70_strong_toxin = max(sigma70_strong_30(:,1)) - max(sigma70_strong_10(:,1)); 

        %% (6) sigma70 weak baseline penalty 
    sigma70_weak_15 = baseline_simulated_data{27}; 
    sigma70_weak_30 = baseline_simulated_data{28}; 
    sigma70_weak_no_toxin = max(sigma70_weak_15(:,1) - max(sigma70_weak_30(:,1))); % We don't want expression to decrease in this case

        %% (7) Positive Crosstalk & (8) Negative Crosstalk (organized into vector) & (9) Promotor Strength
            % Positive crosstalk values should be smaller than 0 except for T7weak
            % emptyT7 case 
            % Negative crosstalk values should all be smaller than 0 
    positive_crosstalk_vec = nan(12,1); % This should contain the maximum crosstalk for 4 promotor strength * 3 combinations = 12  
    negative_crosstalk_vec = nan(12,1); 
    max_strong_crosstalk_ratio = 0; % Keep track of largest crosstalk ratio in T7 strong, T7 weak, and sigma70 strong
    max_weak_crosstalk_ratio = 0; % Keep track of largest crosstalk ratio in sigma70 weak 
    for crosstalk_ratio_row_idx = 1:size(simulated_crosstalk_ratios,1)
        for crosstalk_ratio_col_idx = 1:size(simulated_crosstalk_ratios,2)

            simulated_crosstalk_ratio = simulated_crosstalk_ratios{crosstalk_ratio_row_idx,crosstalk_ratio_col_idx}; 
            positive_crosstalk_vec((crosstalk_ratio_row_idx - 1) * size(simulated_crosstalk_ratios,2) + crosstalk_ratio_col_idx) =...
                min(1 - simulated_crosstalk_ratio(1:3)); % Positive crosstalk should occur in the 3 loweset plasmid concentrations
            negative_crosstalk_vec((crosstalk_ratio_row_idx - 1) * size(simulated_crosstalk_ratios,2) + crosstalk_ratio_col_idx) =...
                min(simulated_crosstalk_ratio(4:7) - 1); % Negative crosstalk should occur in the 4 highest plasmid concentrations

            if ~isequal(crosstalk_ratio_col_idx,size(simulated_crosstalk_ratios,2))
                if max(simulated_crosstalk_ratio(1:3)) > max_strong_crosstalk_ratio % Only account for positive crosstalk in lowest 3 plasmid concentrations 
                    max_strong_crosstalk_ratio = max(simulated_crosstalk_ratio); 
                end
            else
                if max(simulated_crosstalk_ratio(1:3)) > max_weak_crosstalk_ratio
                    max_weak_crosstalk_ratio = max(simulated_crosstalk_ratio); 
                end
            end
        end
    end 
    
    penalty_promotor = max_strong_crosstalk_ratio - max_weak_crosstalk_ratio; 

        %% (10) Other penalties
    other_penalty = calculate_other_penalty(simulated_data,species_name_list);

    %% Organize the penalty terms into a vector form 
    unscaled_penalty_vec = [baseline_dev;crosstalk_ratio_dev;T7_strong_toxin;T7_weak_toxin;...
        sigma70_strong_toxin;sigma70_weak_no_toxin;positive_crosstalk_vec;negative_crosstalk_vec;...
        penalty_promotor;other_penalty;]; 

    penalty_term_labels = {'baseline_dev','crosstalk_ratio_dev','T7_strong_baseline_decrease','T7_weak_baseline_decrease',...
        'sigma70_strong_decrease','sigma70_weak_no_decrease','positive_crosstalk','negative_crosstalk',...
        'promotor_strength','other_penalty'}; 
    penalty_term_length = [1,1,1,1,...
        1,1,length(positive_crosstalk_vec),length(negative_crosstalk_vec),...
        1,1]; 
    % Convert large terms using logs 
        % For each baseline data penalty, convert
    log_T7_strong_toxin = log_conversion(T7_strong_toxin); 
    log_T7_weak_toxin = log_conversion(T7_weak_toxin); 
    log_sigma70_strong_toxin = log_conversion(sigma70_strong_toxin); 
    log_sigma70_weak_no_toxin = log_conversion(sigma70_weak_no_toxin); 
    log_other_penalty = log_conversion(other_penalty); 

    scaled_penalty_vec = [log10(baseline_dev);log10(crosstalk_ratio_dev);log_T7_strong_toxin;log_T7_weak_toxin;...
        log_sigma70_strong_toxin;log_sigma70_weak_no_toxin;positive_crosstalk_vec;negative_crosstalk_vec;...
        penalty_promotor;log_other_penalty;];
    

end