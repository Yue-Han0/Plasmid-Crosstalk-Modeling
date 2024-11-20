function [all_fval,all_estimated_params] = fit_RNAdeg_exp_data(Model,RNA_deg_probObject)

    num_iter = 1200; 

        % Create simFunction 
    model_species = {Model.Species.Name}; 
    species_oi_name = 'RNA utrbroc--no_protein'; 
    species_oi_idx = find(strcmp(species_oi_name,model_species)); 
    track_species = [model_species(1:species_oi_idx - 1),model_species(species_oi_idx + 1:end)];
    simFunction = create_simFun_from_problemObject(RNA_deg_probObject,track_species); 
    dosing_information= create_dosing_info_from_problemObject(RNA_deg_probObject); 
    tStart = 0; 
    tEnd = 2940; 
    options = optimoptions('fmincon','Display','iter');
    
        % Preassign 
    all_fval = nan(num_iter,1);
    all_estimated_params = nan(num_iter,length(RNA_deg_probObject.Estimated)); 
    
        % Run Parameter Estimation 
    param_bounds = [RNA_deg_probObject.Estimated.Bounds]; 
    param_bounds = reshape(param_bounds,[2,length(RNA_deg_probObject.Estimated)]);
    param_bounds = param_bounds'; 
    parfor iter = 1:num_iter
        log_transformed_init_params = nan(1,length(RNA_deg_probObject.Estimated));
        init_params = nan(1,length(RNA_deg_probObject.Estimated));
        for est_param_idx = 1:length(RNA_deg_probObject.Estimated)
            log_transformed_init_params (est_param_idx) =  RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(1) +...
                (RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(2) -  RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(1)) * rand;
            init_params(est_param_idx) = exp(log_transformed_init_params(est_param_idx)); 
        end
        try
            [estimated_params,fval,exitflag,output] = fmincon(@(x)calc_dev(x,simFunction,dosing_information,tStart,tEnd,RNA_deg_probObject.Data),init_params,[],[],[],[],param_bounds(:,1),param_bounds(:,2),[],options);
        
            all_fval(iter) = fval; 
            all_estimated_params(iter,:) = estimated_params; 
        catch
            fprintf('iteration #%d failed',iter)
        end
    end


end


