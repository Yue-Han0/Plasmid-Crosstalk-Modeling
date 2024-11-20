% Write an optimization problem in Gurobi, with the objectives as hard
% constraints and no bounds on parameters; tighten the bounds for
% parameters and see which are important 

% This is a first-pass try to see if the idea works 


clear
clc

currentpath = pwd; 
addpath(genpath(currentpath))
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))

num_iter = 100;

%% A function that takes parameters as input and output penalty terms 

% all_penalty_terms = wrapper_calcQualObj_from_params(params); 
%%  Use sampling results to inform bounds on constraints
sampling_result_file = load('test_save_files/202408_high_res_penalty_valid_terms.mat','valid_high_res_penalty_terms'); 
high_res_penalty_terms = sampling_result_file.valid_high_res_penalty_terms; 

median_baseline_penalty_value = median(high_res_penalty_terms(:,1));
median_crosstalk_ratio_penalty_value = median(high_res_penalty_terms(:,2));
median_other_penalty_value = median(high_res_penalty_terms(:,end)); 

positive_crosstalk_goals = -0.05 .* ones(1,12);
negative_crosstalk_goals = -0.05 .* ones(1,12); 
penalty_goals = [median_baseline_penalty_value,median_crosstalk_ratio_penalty_value,-1e-04,-1e-04,-1e-04,-1e-04,positive_crosstalk_goals,negative_crosstalk_goals,-1e-04,median_other_penalty_value];


%% Can also write in fmincon if gurobi doesn't work 

    % Get information on parameters to be estimated 
group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
    'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
    'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
    'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
param_info_path = 'parameters_unbounded.xlsx';
problemObject = setProblemObject_v3(group_description,[],param_info_path); 
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

% Try various initial points 
all_opt_params = nan(num_iter,num_params);
all_fval = nan(num_iter,1);
all_exitflag = nan(num_iter,1);
all_output = cell(num_iter,1); 

parfor iter = 1:num_iter

    if isequal(iter,1)
        init_params = [problemObject.Estimated.InitialValue];
    else
        log_init_params = 8 * rand(1,num_params);
        init_params = 10 .^ log_init_params;
    end
    try

        [opt_params,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) dumb_obj_func(x,problemObject,dosing_information,simFunction),...
            init_params,[],[],[],[],params_lb,params_ub,@(x) qual_obj_constraint(x,problemObject,dosing_information,penalty_goals),options);
    
        % Assign 
        all_opt_params(iter,:) = opt_params;
        all_fval(iter) = fval;
        all_exitflag(iter) = exitflag;
        all_output{iter} = output; 
    catch
        fprintf('Optimization failed with iter #%d',iter)
    end
end

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


%% Define the gurobi model 

    % This might not work since it's probably not possible to modify the
    % optimization problem constantly using gurobi 
% 
%     % Build model
% model.modelname  = 'test_feasibility';
% model.modelsense = 'min';
% 
%     % Add parameters as variables in the model 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject = setProblemObject_v3(group_description,[],param_info_path); 
% estimated_param_name_list = {problemObject.Estimated.Name};
% num_params = length(estimated_param_name_list); 
% model.lb = zeros(1,num_params);
% model.ub = inf(1,num_params); 
% model.obj = zeros(1,num_params); % Have a null objective function here - just checking feasibility
% model.vtype = 'C';
% 
% model.A = sparse(zeros(32,length(model.lb))); 
% model.sense = repmat('>',32,1); 
% model.rhs = zeros(32,1);  % Temporary rhs values 
% 
% % Set optimization parameters
% params.outputflag = 1;
% 
% % Run the optimization and pass the callback function
% % env = Env("gurobi.log"); 
% result = gurobi(@(model, result) gurobi_callback(model, result), params);
% 
% % Display the results
% disp(result.x);
% 
% % Define the callback function
% function model = gurobi_callback(model, result)
% 
%     if strcmp(result.status, 'LOADED') || strcmp(result.status, 'OPTIMAL')
%         params = result.x; % Get current parameter values
% 
%         % Calculte penalty terms 
%         all_penalty_terms = wrapper_calcQualObj_from_params(params); 
% 
%         % Update rhs 
%         model.rhs = all_penalty_terms;
% 
%     end
% 
% end








