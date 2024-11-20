function [data_matrix_imputed_mean,full_sim_time_course_mean,init_cond_eff_imputed,init_cond_eff_full_sim,full_sim_timeVec] = extrap_mRNA_deg_cont(timeVec_table,data_table,init_cond,plot_flag)
  

    % Extract timeVec
    timeVec = timeVec_table.Time;

    % Add 4 minutes to mimic the actual time 
    time_delay = 240; 
    timeVec_actual = timeVec + time_delay; 

    % For each sample, remove the initially increased part, 
    % assume 1st-order kinetics degradation and extrapolate to time-zero
        % It doesn't seem like individual time-course would work well - use
        % moving average
    table_variables = data_table.Properties.VariableNames;
        % Convert table to matrix, remove index column
    data_matrix = table2array(data_table); 
    data_matrix = data_matrix(:,2:end);
        % Calculate moving average (k = 3) 
    movmean_num = 5; 
    data_matrix_moving_avg = movmean(data_matrix,movmean_num,1); 
        % Calcuate autocorrelation 
    data_matrix_moving_avg_timestepp1 = [data_matrix_moving_avg(2:end,:);...
        data_matrix_moving_avg(end,:) - 0.01 .* ones(1,size(data_matrix_moving_avg,2))];  % add a last row that's guaranteed to be decreasing
    data_autocorrelation = data_matrix_moving_avg_timestepp1 - data_matrix_moving_avg;
        % Find the last timepoint that increases over time and remove any
        % data prior to that timepoint
    data_matrix_filtered = nan(size(data_matrix)); 
    data_timeVec_filtered = nan(size(data_matrix)); 
    for col_idx = 1:size(data_matrix_filtered,2)
        start_remove_idx_list = find(data_autocorrelation(:,col_idx) > 0);
        max_row_keep_idx = max(start_remove_idx_list);
        if isempty(max_row_keep_idx) % if monotonically decreasing, add the entire time course
            data_matrix_filtered(:,col_idx) = data_matrix(:,col_idx); 
            data_timeVec_filtered(:,col_idx) = timeVec_actual; 
        else 
            data_matrix_filtered(max_row_keep_idx+1:end,col_idx) = data_matrix(max_row_keep_idx+1:end,col_idx); 
            data_timeVec_filtered(max_row_keep_idx+1:end,col_idx) = timeVec_actual(max_row_keep_idx+1:end); 
        end
    end

    % Fit the filtered matrix with an exponential decay function
    % First-order decay: ln(ct) = ln(c0) - kt
        % Recored fitted c0, k, and R^2
    C0_vec = nan(1,size(data_matrix,2));
    k_vec = nan(1,size(data_matrix,2)); 
    R_squared_vec =  nan(1,size(data_matrix,2)); 
    for col_idx2 = 1:size(data_matrix_filtered,2)
        data_matrix_col_oi = data_matrix_filtered(:,col_idx2); 
        data_matrix_col_oi_eff = data_matrix_col_oi(~isnan(data_matrix_col_oi)); 
        if length(data_matrix_col_oi_eff) > 10 % Ignore time-courses with fewer than 10 timepoints 
            lnCt = log(data_matrix_col_oi_eff); 
            t = data_timeVec_filtered(~isnan(data_timeVec_filtered(:,col_idx2)),col_idx2); 
            [lr_coeff,S] = polyfit(t,lnCt,1); 
            R_squared = 1 - (S.normr/norm(lnCt - mean(lnCt)))^2;
            C0_vec(col_idx2) = exp(lr_coeff(2));
            k_vec(col_idx2) = -lr_coeff(1);
            R_squared_vec(col_idx2) = R_squared; 
        end
    
    end

    % If LR is reasonable, fill the missing data with imputed values 
    data_matrix_imputed = data_matrix_filtered; 
    full_sim_timeVec = 0:30:timeVec_actual(end); 
    full_sim_time_course = nan(length(full_sim_timeVec),size(data_matrix_filtered,2)); 
    for col_idx3 = 1:size(data_matrix_filtered,2)
        if ~isnan(R_squared_vec(col_idx3)) && any(isnan(data_matrix_filtered(:,col_idx3)))
            % Calculate fitted time-course
            sim_time_course = exp(log(C0_vec(col_idx3)) - k_vec(col_idx3) .* timeVec_actual);
    
            % Fill in missing value
            data_matrix_imputed(isnan(data_matrix_imputed(:,col_idx3)),col_idx3) = sim_time_course(isnan(data_matrix_imputed(:,col_idx3)));
    
        end
        if ~isnan(R_squared_vec(col_idx3))
            % Calculated fully simulated time-course
            full_sim_time_course(:,col_idx3) =  exp(log(C0_vec(col_idx3)) - k_vec(col_idx3) .* full_sim_timeVec);
        end
    end

    if plot_flag
        % Compare imputed value & simulated value with original values 
        figure;
        for col_idx4 = 1:3:size(data_matrix_imputed,2)
            subplot(3,3,ceil(col_idx4 / 3))
            for rep = 1:3
                %     % Plot initial actual time-course
                % plot(timeVec_actual(1:10),data_matrix(1:10,col_idx4 + rep - 1),'LineWidth',1.5,'Color','k')
                % hold on 
                % plot(timeVec_actual(1:10),data_matrix_imputed(1:10,col_idx4 + rep - 1),'LineWidth',1.5,'LineStyle','--','Color','r')
                % plot(full_sim_timeVec(1:20),full_sim_time_course(1:20,col_idx4 + rep - 1),'LineWidth',1.5,'LineStyle','-.','Color','b')
                    % Plot all time-course
                plot(timeVec_actual,data_matrix(:,col_idx4 + rep - 1),'LineWidth',1.5,'Color','k')
                hold on 
                plot(timeVec_actual,data_matrix_imputed(:,col_idx4 + rep - 1),'LineWidth',1.5,'LineStyle','--','Color','r')
                plot(full_sim_timeVec,full_sim_time_course(:,col_idx4 + rep - 1),'LineWidth',1.5,'LineStyle','-.','Color','b')
            end
            title(sprintf('[mRNA]_{init} = %d nM',init_cond(ceil(col_idx4 / 3))))
        end
        legend('Experimental','Imputed','Simulated')
        % sgtitle('')
    end

    % if any of the 3 replicates are present, calculate the average 
    data_matrix_imputed_mean = nan(size(data_matrix_imputed,1),size(data_matrix_imputed,2)/3);
    full_sim_time_course_mean = nan(size(full_sim_time_course,1),size(full_sim_time_course,2)/3);
    init_cond_eff_imputed = []; 
    init_cond_eff_full_sim = []; 
    for col_idx5 = 1:size(data_matrix_imputed,2)/3
            % only calculate the mean if at least one third of data is
            % present 
        if sum(sum(isnan(data_matrix_imputed(:,(col_idx5 - 1) * 3 + 1:col_idx5 * 3)))) <= ...
                0.67 * numel(data_matrix_imputed(:,(col_idx5 - 1) * 3 + 1:col_idx5 * 3))
            data_matrix_imputed_mean(:,col_idx5) = mean(data_matrix_imputed(:,(col_idx5 - 1) * 3 + 1:col_idx5 * 3),2,'omitmissing');
            init_cond_eff_imputed = [init_cond_eff_imputed init_cond(col_idx5)];
        end
        if sum(sum(isnan(full_sim_time_course(:,(col_idx5 - 1) * 3 + 1:col_idx5 * 3)))) <= ...
                0.67 * numel(full_sim_time_course(:,(col_idx5 - 1) * 3 + 1:col_idx5 * 3))
            full_sim_time_course_mean(:,col_idx5) = mean(full_sim_time_course(:,(col_idx5 - 1) * 3 + 1:col_idx5 * 3),2,'omitmissing');
            init_cond_eff_full_sim = [init_cond_eff_full_sim init_cond(col_idx5)];
        end
    end

 


end