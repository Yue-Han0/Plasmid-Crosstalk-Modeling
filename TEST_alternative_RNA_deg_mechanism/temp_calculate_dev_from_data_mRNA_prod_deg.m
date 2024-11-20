function varargout = temp_calculate_dev_from_data_mRNA_prod_deg(paramVec)

    file_struct = load('param_est_run_save/20231019_mRNA_prod_deg_model_fitting_analysis.mat');
    problemObject = file_struct.problemObject;

    sim_function = create_simFun_from_problemObject(problemObject);
    dosing_information = create_dosing_info_from_problemObject(problemObject);
    tStart = 0; 
    tEnd = 21600; 

    [sim_time,sim_data] = sim_function(paramVec',tEnd,dosing_information,tStart:tEnd); 

    % At some point, double-check that the data and dosing information are
    % in the same order 
    [exp_time,exp_data] = process_crosstalk_data_from_source(problemObject.Data,'grouped_data'); 

    % calculate the deviation between experimental vs. simulated data
    dev = 0; 
    for dataset_idx = 1:length(exp_time)
        sim_time_single = sim_time{dataset_idx};
        sim_data_single = sim_data{dataset_idx};

        exp_time_single = exp_time{dataset_idx};
        exp_data_single = exp_data{dataset_idx}; 

        % process sim_time to avoid error in interp1 
        [unique_sim_data_single,unique_idx] = unique(sim_data_single); 
        if ~isequal(length(unique_sim_data_single),length(sim_data_single))
            unique_sim_time_single = sim_time_single(unique_idx); 
            extrap_sim_data_single = interp1(unique_sim_time_single,unique_sim_data_single,exp_time_single,'linear','extrap'); 
        else
            extrap_sim_data_single = interp1(sim_time_single,sim_data_single,exp_time_single,'linear','extrap'); 
        end

        dev = dev + sum((exp_data_single - extrap_sim_data_single).^2);
    end
    if (nargout == 1)
        varargout{1} = dev;
    end

end



