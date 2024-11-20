function [dev,penalty_transition,penalty_promotor] = calculate_crosstalk_ratio_dev(experimental_crosstalk_ratios,simulated_crosstalk_ratios)


    % Calculate the deviation of crosstalk ratios 
    dev = 0; 
    for crosstalk_ratio_row_idx = 1:size(experimental_crosstalk_ratios,1)
        for crosstalk_ratio_col_idx = 1:size(experimental_crosstalk_ratios,2)
            experimental_crosstalk_ratio = experimental_crosstalk_ratios{crosstalk_ratio_row_idx,crosstalk_ratio_col_idx};
            simulated_crosstalk_ratio = simulated_crosstalk_ratios{crosstalk_ratio_row_idx,crosstalk_ratio_col_idx}; 
            dev = dev + sum((simulated_crosstalk_ratio - experimental_crosstalk_ratio).^2); 
        end
    end

    % Calcualte penalty for (1) positive/negative crosstalk ratios (2) weak
    % promotors 
    penalty_transition = 0; 
    max_strong_crosstalk_ratio = 0; % Keep track of largest crosstalk ratio in T7 strong, T7 weak, and sigma70 strong
    max_weak_crosstalk_ratio = 0; % Keep track of largest crosstalk ratio in sigma70 weak 
    for crosstalk_ratio_row_idx = 1:size(simulated_crosstalk_ratios,1)
        for crosstalk_ratio_col_idx = 1:size(simulated_crosstalk_ratios,2)
            simulated_crosstalk_ratio = simulated_crosstalk_ratios{crosstalk_ratio_row_idx,crosstalk_ratio_col_idx}; 
            penalty_transition = penalty_transition + min(1 - simulated_crosstalk_ratio(1:3)) + min(simulated_crosstalk_ratio(4:7) -1); 
            if ~isequal(crosstalk_ratio_col_idx,size(simulated_crosstalk_ratios,2))
                if max(simulated_crosstalk_ratio) > max_strong_crosstalk_ratio
                    max_strong_crosstalk_ratio = max(simulated_crosstalk_ratio(1:3)); 
                end
            else
                if max(simulated_crosstalk_ratio) > max_weak_crosstalk_ratio
                    max_weak_crosstalk_ratio = max(simulated_crosstalk_ratio(1:3)); 
                end
            end
        end
    end

    penalty_promotor = max_strong_crosstalk_ratio - max_weak_crosstalk_ratio; 
end