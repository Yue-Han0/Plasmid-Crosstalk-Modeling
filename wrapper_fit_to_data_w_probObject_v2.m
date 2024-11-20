function wrapper_fit_to_data_w_probObject_v2(problemObject,save_path_id,override_result_path)
    % Perform parameter estimation with plasmid crosstalk data 
    % Update from v2: sample in the log space 
    
    % Assumes data is loaded from a table 

    % save_path_id: a string for experiment description 
    % problemObject: Simbiology problem object

%%%%%%%%%setting ode solver to ode15s%%%%%%%%%%%
% model_solver = getconfigset(problemObject.Model);
% model_solver.SolverType = "ode15s";
% problemObject.Options = optimoptions('lsqnonlin','Display','iter'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define save path and pre-assigned values 
numRuns = 100;
if ~exist("override_result_path",'var')
    time = datetime; 
    save_path_stem = sprintf('param_est_run_save/2024%02d%02d_param_est_run%02d%02d_%s',time.Month,time.Day,time.Hour,time.Minute,save_path_id); 
else
    save_path_stem = override_result_path; 
end

problemObject.Options = optimoptions('lsqnonlin','Display','iter'); 

% Run fit 
if ~exist("override_result_path",'var')
    all_fitResults = cell(numRuns,1); 
    start_idx = 1;
else
    result_file = load(override_result_path);
    all_fitResults = result_file.all_fitResults; 
    num_completed_run = nnz(~cellfun(@isempty,all_fitResults));
    start_idx = num_completed_run + 1; 
end
parfor iter = start_idx:numRuns
    % Modify initial conditions, sample within bounds
    success_flag = 0; 
    num_trial = 0; 
    while ~success_flag % Adding this to avoid parameter estimation failure due to bad initial point 
        rng(iter + num_trial * numRuns)
        temp_problemObject = problemObject; 
        % Use initial parameter values for first run 
        if ~isequal(iter,1)
            % Then sample initial parameter values within parameter bounds
            % in the log space 
            for est_param_idx = 1:length(temp_problemObject.Estimated)
                temp_problemObject.Estimated(est_param_idx).InitialTransformedValue =  temp_problemObject.Estimated(est_param_idx).TransformedBounds(1) + ...
                    (temp_problemObject.Estimated(est_param_idx).TransformedBounds(2) -  temp_problemObject.Estimated(est_param_idx).TransformedBounds(1)) * rand;
                temp_problemObject.Estimated(est_param_idx).InitialValue =  exp(temp_problemObject.Estimated(est_param_idx).InitialTransformedValue);
            end
        end
        try
            [fitResults,simdataI] = fit(temp_problemObject); % A fit problem object 
            success_flag = 1;
        catch
            success_flag = 0; 
            num_trial = num_trial + 1; 
        end
    end
    
    all_fitResults{iter,1} = fitResults; 
    % save(sprintf('%s.mat',save_path_stem),'problemObject','all_fitResults','-v7.3'); % Save along the way
end

% Save results 
try
    save(sprintf('%s.mat',save_path_stem),'problemObject','all_fitResults','-v7.3');
catch
    for i = 1:numRuns
        fitResults = all_fitResults{i,1}; 
        save(sprintf('%s_fitResult%d.mat',save_path_stem,i),'fitResults','problemObject','-v7.3');
    end
end



end