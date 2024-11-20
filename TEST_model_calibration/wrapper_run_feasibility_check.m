function wrapper_run_feasibility_check(mechanism,iter)


    currentpath = pwd; 
    addpath(genpath(currentpath))
    rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))


    %% Load penalty goal file
    penalty_goal_file = load('test_save_files/high_resolution_penalty_goals.mat');
    penalty_goals = penalty_goal_file.penalty_goals; 

    %% Can also write in fmincon if gurobi doesn't work 
    
        % Get information on parameters to be estimated 
    group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
        'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
        'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
        'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
    param_info_path = 'parameters_v2_newMech.xlsx';

    if strcmp(mechanism,'TestBindingSites')
        problemObject = setProblemObject_v3_TestBindingSite(group_description,[],param_info_path);
    elseif strcmp(mechanism,'TestPrototype')
        problemObject = setProblemObject_v3_TestPrototype(group_description,[],param_info_path); 
    else
        error('Unspecified Testing Mechanism')
    end

    dosing_information = create_dosing_info_from_problemObject(problemObject); 
    estimated_param_name_list = {problemObject.Estimated.Name};
    num_params = length(estimated_param_name_list); 
    
        % Create simFunction
            % Reorganize observables such that sfGFP comes first 
    model = problemObject.Model; 
    species_name_list = {model.Species.Name}; 
    remove_idx = strcmp('protein sfGFP*',species_name_list); 
    additional_track_species_list = species_name_list(~remove_idx); 
    observables = [{'protein sfGFP*'},additional_track_species_list];
    
    % % Calculate time-course data using the current set of parameters 
    % Alternatively 
    simFunction = create_simFun_from_problemObject(problemObject,observables(2:end)); 
    
    params_lb = zeros(num_params,1);
    params_ub = inf(num_params,1);
    options = optimoptions('fmincon','Display','iter',"EnableFeasibilityMode",true,...
        "SubproblemAlgorithm","cg");
    
    num_trial = 1;
    successFlag = false; 
    while ~successFlag && num_trial < 50
        if isequal(iter,1)
            init_params = [problemObject.Estimated.InitialValue];
        else
            rng(iter * 50 + num_trial)
            log_init_params = -8 + 16 * rand(1,num_params);
            init_params = 10 .^ log_init_params;
        end
        try
            [opt_params,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) dumb_obj_func(x,problemObject,dosing_information,simFunction),...
                init_params,[],[],[],[],params_lb,params_ub,@(x) qual_obj_constraint(x,problemObject,dosing_information,penalty_goals),options);
            successFlag = true; 
        catch
            num_trial = num_trial + 1; 
            
        end
    end
    save(sprintf('param_est_run_save/20241010_%s_Feasibility_Check_iter%d_updated.mat',mechanism,iter));
 

    % A function that construct constraints from calculated penalty terms 
    function [c,ceq] = qual_obj_constraint(params,problemObject,dosing_information,penalty_goals)
    
        all_penalty_terms = wrapper_calcQualObj_from_params(params,problemObject,dosing_information); 
    
        c = all_penalty_terms - penalty_goals; 
        ceq = []; 
    
    end

    % A dumb function for objectives (minimizing SSE between experimental and
    % simulated data) 
    function obj = dumb_obj_func(params,problemObject,dosing_information,simFunction)
        tStart = 0; 
        tEnd = 21600; % Modify this later to draw info from data  
    
        experimental_grouped_data = problemObject.Data; 
        [simulated_time,simulated_data] = simFunction(params,tEnd,dosing_information,tStart:tEnd);
        [experimental_time,experimental_data] = process_crosstalk_data_from_source(experimental_grouped_data,'grouped_data'); 
    
        % Extrapolate simulated data to timeVec of the experimental data
        obj = 0; 
        for dataset_idx = 1:length(simulated_data)
    
            experimental_time_single = experimental_time{dataset_idx}; 
            experimental_data_single = experimental_data{dataset_idx}; 
            simulated_time_single = simulated_time{dataset_idx}; 
            simulated_data_single = simulated_data{dataset_idx}; 
    
            extrap_simulated_data_single = interp1(simulated_time_single(1:10:end),simulated_data_single(1:10:end,1),...
                experimental_time_single,'linear','extrap'); 
            obj = obj + sum((experimental_data_single - extrap_simulated_data_single(:,1)).^2); 
        end
    end


end