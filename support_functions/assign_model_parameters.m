function modified_Mobj = assign_model_parameters(Mobj,param_names,param_vals)
    % Assign parameter values to parameter names in a Simbiology model

    % Mobj - A Simbiology object
    % param_names - A cell array with parameter names aligning those in
    % Mobj.Parameters
    % param_vals - A vector with parameter values 

    modified_Mobj = Mobj; 
    for param_idx = 1:length(param_names)
        param_name = param_names{param_idx};
        param_val = param_vals(param_idx); 
        param_obj_oi = sbioselect(modified_Mobj,'Type','Parameter','Name',param_name); 
        if ~isempty(param_obj_oi)
            set(param_obj_oi,'Value',param_val)
        end
    end


end