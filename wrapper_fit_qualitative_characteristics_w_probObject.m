function wrapper_fit_qualitative_characteristics_w_probObject(problemObject,sampled_params,save_path_id)

    %% Dumb variables 
    conc_vec = [0.5,1,2.5,5,10,15,30]; 
    mode = 'PE';
    numRuns = size(sampled_params,1); 
    num_conc = length(conc_vec);
    num_promotor = 4;
    obj_string = 'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + crosstalk_ratio_transition_penalty + crosstalk_ratio_promotor_strength_penalty + log10(other_penalty)';
    
    % Set problem object and get experimental data 
    experimental_grouped_data = problemObject.Data;
    Mobj = problemObject.Model; 
    param_bounds = [problemObject.Estimated.Bounds]; 
    param_bounds = reshape(param_bounds,[2,length(problemObject.Estimated)]);
    param_bounds = param_bounds'; 
    lb = param_bounds(:,1);
    ub = param_bounds(:,2); 
    
    %% Create simFunction and track all species 
    GFP_species_idx = find(strcmp({Mobj.Species.Name},'protein sfGFP*'));
    all_species_names = {Mobj.Species.Name}; 
    additional_track_species = [all_species_names{1:GFP_species_idx-1},all_species_names(GFP_species_idx+1:end)]; 
    species_name_list = [{'GFP'},additional_track_species];
    simFunction = create_simFun_from_problemObject(problemObject,additional_track_species);
    
    options = optimoptions('fmincon','Display','iter');
    % rxn_oi_list = {'[CUTP:AGTP:RNAP:DNA pJ23105--utrGFP--sfGFP] -> [term_RNAP:DNA pJ23105--utrGFP--sfGFP] + [RNA utrGFP--sfGFP]'}; 
    % param_oi_list = {'TXTL_PJ23105_RNAPbound_R'}; 
    % modify_value_magnitude_list = -4:1:4; 
    % modify_values = 10.^ modify_value_magnitude_list;
    % modify_value_list = {modify_values}; 
    % sensitivity_goal_list = {1}; 
    % sensitivity_options = struct('rxn_oi_list',rxn_oi_list,'param_oi_list',param_oi_list,'modify_value_list',modify_value_list,...
    %     'sensitivity_goal_list',sensitivity_goal_list); 
    sensitivity_options = [];
    
    %% Run parameter estimation 
        % Add in continue-run mechanism
    if ~exist(save_path_id,'file')
        all_fitResults = cell(numRuns,1); 
        penalty_term_label = {'baseline_data_dev','crosstalk_ratio_dev','log_baseline_data_penalty',...
            'crosstalk_ratio_transition_penalty','crosstalk_ratio_promotor_strength_penalty','promotor_sensitivity_penalty','other_penalty'}; 
        all_penalty_term_init = nan(numRuns,length(penalty_term_label)); 
        all_penalty_term_final = nan(numRuns,length(penalty_term_label)); 
        start_idx = 1; 
    else
        result_file = load(save_path_id);
        all_fitResults = result_file.all_fitResults; 
        penalty_term_label = result_file.penalty_term_label;
        all_penalty_term_init = result_file.all_penalty_term_init; 
        all_penalty_term_final = result_file.all_penalty_term_final; 
        num_completed_run = nnz(~cellfun(@isempty,all_fitResults));
        start_idx = num_completed_run + 1; 
    end

    for iter = start_idx:numRuns

        sampled_params_single = sampled_params(iter,:); 
        % Initial penalty term values 
        [~,baseline_data_dev_init,crosstalk_ratio_dev_init,log_baseline_data_penalty_init,crosstalk_ratio_transition_penalty_init,crosstalk_ratio_promotor_strength_penalty_init,promotor_sensitivity_penalty_init,other_penalty_init] = ...
            wrapper_calculate_obj(sampled_params_single,obj_string,simFunction,problemObject,experimental_grouped_data,species_name_list,sensitivity_options,mode,num_conc,num_promotor); 
        % [estimated_params,resnorm,residual,exitflag,output] = lsqnonlin(@(x) wrapper_calculate_obj(x,obj_string,simFunction,problemObject,experimental_grouped_data,species_name_list,mode,num_conc,num_promotor),init_params,lb,ub,options);
        [estimated_params,fval,exitflag,output] = fmincon(@(x) wrapper_calculate_obj(x,obj_string,simFunction,problemObject,experimental_grouped_data,species_name_list,sensitivity_options,mode,num_conc,num_promotor),sampled_params_single,[],[],[],[],lb,ub,[],options);
        fitResult.estimated_params = estimated_params;

        % Final penalty term values 
        [~,baseline_data_dev_final,crosstalk_ratio_dev_final,log_baseline_data_penalty_final,crosstalk_ratio_transition_penalty_final,crosstalk_ratio_promotor_strength_penalty_final,promotor_sensitivity_penalty_final,other_penalty_final] = ...
            wrapper_calculate_obj(estimated_params,obj_string,simFunction,problemObject,experimental_grouped_data,species_name_list,sensitivity_options,mode,num_conc,num_promotor); 
        fitResult.fval = fval; 
        fitResult.exitflag = exitflag;
        fitResult.output = output; 
        all_fitResults{iter,1} = fitResult; 
        all_penalty_term_init(iter,:) = [baseline_data_dev_init,crosstalk_ratio_dev_init,log_baseline_data_penalty_init,crosstalk_ratio_transition_penalty_init,crosstalk_ratio_promotor_strength_penalty_init,promotor_sensitivity_penalty_init,other_penalty_init]; 
        all_penalty_term_final(iter,:) = [baseline_data_dev_final,crosstalk_ratio_dev_final,log_baseline_data_penalty_final,crosstalk_ratio_transition_penalty_final,crosstalk_ratio_promotor_strength_penalty_final,promotor_sensitivity_penalty_final,other_penalty_final]; 

        % Save results 
        save(save_path_id,'simFunction','all_fitResults','problemObject','penalty_term_label','all_penalty_term_final','all_penalty_term_init');

    end
    
    


end