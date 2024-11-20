function wrapper_analyze_fit_result_v2(all_fitResults,problemObject,fit_result_path,metric_name,plot_fitting,plot_param_dist) 
    
    
    numRuns = length(all_fitResults); 

    % Get the defined metric and estimated parameter summary 
    metric_summary = nan(numRuns,1);
    
    for iter = 1:numRuns
        fitResult = all_fitResults{iter,1};  
        eval(sprintf('metric_summary(iter,1) = fitResult.%s;',metric_name))
        est_param_table = fitResult.ParameterEstimates; 
        if ~exist('estimated_param_summary','var')
            estimated_param_summary = est_param_table;
            estimated_param_summary = removevars(estimated_param_summary,'StandardError');
            estimated_param_summary = renamevars(estimated_param_summary,'Estimate','Estimate_1');
        else
            new_col_name = sprintf('Estimate_%d',iter); 
            estimated_param_summary.(new_col_name) = est_param_table.Estimate; 
        end
    end
    
    % Save this summary into the result structure
    try 
        save(fit_result_path,'all_fitResults','problemObject','estimated_param_summary','metric_summary','metric_name')
    catch
        save(fit_result_path,'problemObject','estimated_param_summary','metric_summary','metric_name')
    end
    
    % Get the lower and upper bound from problemObject 
    estimated_param_obj_all = problemObject.Estimated; 
    parameter_bounds = nan(length(estimated_param_obj_all),2); 
    for est_idx = 1:length(estimated_param_obj_all)
        estimated_param_obj = estimated_param_obj_all(est_idx); 
        parameter_bounds(est_idx,:) = estimated_param_obj.Bounds; 
    end
    

    % Plot the result with highest log likelihood  
    if exist('plot_fitting','var') && plot_fitting
        if strcmp(metric_name,'LogLikelihood')
            [~,opt_idx] = max(metric_summary); 
        else
            [~,opt_idx] = min(metric_summary); 
        end
        opt_fitResult = all_fitResults{opt_idx,1}; 
        plot(opt_fitResult)
        saveas(gcf,sprintf('%s_fitting',fit_result_path),'png'); 
    end
    
    % Plot the distribution of parameters in numRuns
    if exist('plot_param_dist','var') && plot_param_dist
        figure;
        sgtitle('Distribution of estimated parameters');
        num_params = size(estimated_param_summary,1); 
        for param_idx = 1:num_params
            subplot(ceil(sqrt(num_params)),ceil(sqrt(num_params)),param_idx)
            param_name = estimated_param_summary{param_idx,'Name'}{1,1}; 
            param_name = strrep(param_name,'_',' ');
            param_val = nan(numRuns,1); 
            for iter = 2:numRuns
                param_val(iter,1) = table2array(estimated_param_summary(param_idx,sprintf('Estimate_%d',iter))); 
            end
            histogram(param_val,30)
            hold on 
            xlim(parameter_bounds(param_idx,:)); 
            title(param_name)
        end
        saveas(gcf,sprintf('%s_param_dist',fit_result_path),'png'); 
    end
    
    





end