function [rxn_name_list,rxn_time_course] = get_reaction_time_course(model,estimated_parameter_name_list,estimated_parameter_value_list,simulated_time,simulated_data,simulated_dataNames)
    % Calculate time course for model reactions based on species
    % time-course and kinetic rate laws 

    %% Input 
    % model: a Simbiology model object with reaction and parameter objects 
    % fitResult: a Simbiology OptimResults object with estimated parameters
    % Data: a SimData object with Time, Data, and DataNames

    %% Output 
    % rxn_name_list: #reactions * 1 cell with reaction stoichiometry
    % rxn_time_course: #t timepoints * # reactions with reaction rates 
    
    %% Get model reactions, parameters, and fitted parameters 
    model_reaction_obj_list = model.Reactions; 
    model_parameter_obj_list = model.Parameters; 
    
    %% Assign time-course to each species 
    
    for species_idx = 1:length(simulated_dataNames)
        species_name = simulated_dataNames{species_idx};
            % Replace some characters to avoid assigning error
        if strcmp(species_name,'protein sfGFP*')
            species_name = 'protein_sfGFP_m';
        end
        species_name = strrep(species_name,' ','');
        species_name = strrep(species_name,'--',''); 
        species_name = strrep(species_name,':',''); 
    
        eval(sprintf('%s = simulated_data(:,species_idx);',species_name)); 
    
    end
    
    %% Assign parameter values 
        % Existing parameters in problemObject
    for param_idx_1 = 1:length(model_parameter_obj_list)
        parameter_obj = model_parameter_obj_list(param_idx_1); 
        parameter_name = parameter_obj.Name; 
        eval(sprintf('%s = parameter_obj.Value;',parameter_name)); 
    end
        % Update these based on fitResult
    for param_idx_2 = 1:length(estimated_parameter_name_list)
        parameter_name = estimated_parameter_name_list{param_idx_2}; 
        eval(sprintf('%s = estimated_parameter_value_list(param_idx_2);',parameter_name)); 
    end
    
    %% Assign a matrix for reaction rate calculation 
    rxn_time_course = nan(length(simulated_time),length(model_reaction_obj_list)); 
    
    for rxn_idx = 1:length(model_reaction_obj_list)
        reaction_obj = model_reaction_obj_list(rxn_idx); 
        rxn_rate_string = reaction_obj.ReactionRate; 
            % Replace strings accordingly compared to previous assignments 
        rxn_rate_string = strrep(rxn_rate_string,' ','');
        rxn_rate_string = strrep(rxn_rate_string,'--','');
        rxn_rate_string = strrep(rxn_rate_string,':','');
        rxn_rate_string = strrep(rxn_rate_string,'protein_sfGFP*','protein_sfGFP_m');
        rxn_rate_string = strrep(rxn_rate_string,'[','');
        rxn_rate_string = strrep(rxn_rate_string,']','');
            % Get .* & ./ for correct dimension
        rxn_rate_string = strrep(rxn_rate_string,'*','.*'); 
        rxn_rate_string = strrep(rxn_rate_string,'/','./'); 

        rxn_time_course(:,rxn_idx) = eval(rxn_rate_string); 
    
    
    end
    rxn_name_list = {model_reaction_obj_list.Reaction}'; 
end