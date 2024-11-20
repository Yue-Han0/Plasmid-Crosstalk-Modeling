function Model = modify_Model_RNAdeg_mechanism(Model,kinetics)

    switch kinetics

        case 'TestPrototype'

            % Let's first add new species
            RNase_bound_functional_broc = addspecies(Model,'RNase_bound_f_broc');
            RNase_bound_functional_GFP = addspecies(Model,'RNase_bound_f_GFP');
            RNase_bound_deactivated = addspecies(Model,'RNase_bound_d');
            mRNA_deact = addspecies(Model,'mRNA_deact');

            % Collect all mRNA species in the model
            all_species_name_list = {Model.Species.Name};

            % Identify all RNA species present 
            RNA_species_idx_list = contains(all_species_name_list,'RNA utr') & ~contains(all_species_name_list,':'); 
            RNA_species_name_list = all_species_name_list(RNA_species_idx_list); 

            % For all RNase bidning reactions except with reporter, change the product into
            % RNase_bound_d; for reporter RNase binding reaction, change
            % the product into RNase_bound_f
            deg_rxns_to_be_deleted = {}; 
            deg_rxns_to_be_added = {}; 
            for rxn_idx = 1:length(Model.Reactions)

                % Get both sides on the chemical equation 
                reaction_selected = Model.Reactions(rxn_idx); 
                chemical_reaction = reaction_selected.Reaction; 
                divider_idx = strfind(chemical_reaction,'->'); 
                rxn_lhs = chemical_reaction(1:divider_idx); 
                rxn_rhs = chemical_reaction(divider_idx + 1:end); 

                % for RNase binding reactions, track them and replace the
                % product by RNase_bound_f
                if contains(rxn_lhs,' RNase ') && contains(rxn_rhs,':RNase')  

                    % % Add to the tracking pool
                    % binding_rxns_to_be_added{end + 1} = reaction_selected;

                        % Get product name and remove product 
                    product_name = get(reaction_selected.product); 
                    rmproduct(reaction_selected, product_name.Name); 

                        % Replace with RNase_bound
                    if contains(rxn_lhs,'broc')
                        addproduct(reaction_selected,RNase_bound_functional_broc,1); 
                    elseif contains(rxn_lhs,'sfGFP')
                        addproduct(reaction_selected,RNase_bound_functional_GFP,1); 
                    else
                        addproduct(reaction_selected,RNase_bound_deactivated,1);
                    end

                end

                % Find RNA degradation reaction to add to tracking for deletion later
                if contains(rxn_lhs,':RNase') && contains(rxn_rhs,' RNase ')
                    if ~contains(rxn_lhs,'broc') && ~contains(rxn_lhs,'sfGFP')
                        deg_rxns_to_be_deleted{end + 1} = reaction_selected;
                    else
                        deg_rxns_to_be_added{end + 1} = reaction_selected; 

                        % Get reactant name and remove reactant 
                        reactant_name = get(reaction_selected.Reactants);
                        rmreactant(reaction_selected,reactant_name.Name); 
    
                        % Replace with RNase_bound
                        if contains(rxn_lhs,'broc')
                            addreactant(reaction_selected,RNase_bound_functional_broc,1); 
                        elseif contains(rxn_lhs,'sfGFP')
                            addreactant(reaction_selected,RNase_bound_functional_GFP,1); 
                        end
    
                        % Add 2 deactivated mRNA as product 
                        addproduct(reaction_selected,mRNA_deact,2)
    
                        % Remove AGMP/CUMP 
                        rmproduct(reaction_selected,'AGMP');
                        rmproduct(reaction_selected,'CUMP');
                    end

                end

                % For reactions where mRNA is produced, add a null ->
                % mRNA_deact reaction with the same reaction rate 
                all_rxn_rate_string = {}; 
                all_rxn_rate_string_species_ref = {}; 
                for RNA_species_idx = 1:length(RNA_species_name_list)
                    RNA_species_name = RNA_species_name_list{RNA_species_idx};
                    RNA_species_name_for_model = strrep(RNA_species_name,' ','_');
                    RNA_species_name_for_model = strrep(RNA_species_name_for_model,'-','');


                    if contains(rxn_rhs,strcat('[',RNA_species_name,']')) % then it's a production reaction for species 
                        selected_rxn_rate_string = reaction_selected.ReactionRate; 
                        rxn_string_var_name = sprintf('%s_prod_rxn_rate_string',RNA_species_name_for_model); 
                        if ~exist(rxn_string_var_name,'var')
                            eval(sprintf('%s = selected_rxn_rate_string;',rxn_string_var_name)); 
                            all_rxn_rate_string{end + 1} = rxn_string_var_name;
                            all_rxn_rate_string_species_ref{end + 1} = RNA_species_name_for_model; 
                        else
                            eval(sprintf('%s = strcat(%s,''+'',selected_rxn_rate_string);',rxn_string_var_name,rxn_string_var_name)); 
                        end

                    end
                end

            end

            for prod_rxn_idx = 1:length(all_rxn_rate_string_species_ref)

                prod_species = all_rxn_rate_string_species_ref{prod_rxn_idx};
                prod_rxn_rate_string = all_rxn_rate_string{prod_rxn_idx}; 

                eval(sprintf('%s_prod_rxn = addreaction(Model,''null -> mRNA_deact'');',prod_species));
                eval(sprintf('%_prod_rxn_kineticLaw = addkineticlaw(%s_prod_rxn,''Unknown'');',prod_species,prod_species)); 
                eval(sprintf('set(%s_prod_rxn,''ReactionRate'',prod_rxn_rate_string);',prod_species)); 

            end
    
            % Delete binding reactions with non-reporter mRNA species to
            % avoid double-counting 
            for rxn_idx = 1:length(Model.Reactions)

                % Get both sides on the chemical equation 
                reaction_selected = Model.Reactions(rxn_idx); 
                chemical_reaction = reaction_selected.Reaction; 
                divider_idx = strfind(chemical_reaction,'->'); 
                rxn_lhs = chemical_reaction(1:divider_idx); 
                rxn_rhs = chemical_reaction(divider_idx + 1:end); 

                if contains(rxn_rhs,'RNase_bound_d') && ~contains(rxn_lhs,'mRNA_deact')
                    deg_rxns_to_be_deleted{end + 1} = reaction_selected; 
                end
            end

            % Delete the reactions 
            for delete_rxn_idx = 1:length(deg_rxns_to_be_deleted)
                delete_rxn_selected = deg_rxns_to_be_deleted{delete_rxn_idx}; 
                delete(delete_rxn_selected)
            end

                % Add parameters to represent degradation kinetic constant
                % for deactivated mRNAs; these parameters are technically
                % time-dependent but the average will be taken 
            TXTL_RNAdeg_kc_2 = addparameter(Model,'TXTL_RNAdeg_kc_2',1); 
            TXTL_RNAdeg_kc_1 = addparameter(Model,'TXTL_RNAdeg_kc_1',1); 
            TXTL_RNAdeg_kc_0 = addparameter(Model,'TXTL_RNAdeg_kc_0',1); 
                        % Add degradation parameters to the model
            TXTL_RNAdeg_F = addparameter(Model,'TXTL_RNAdeg_F',1); 
            TXTL_RNAdeg_R = addparameter(Model,'TXTL_RNAdeg_R',1); 

            
              % Add binding reaction with mRNA_deactivated 
            RNase_binding_mRNAdeact = addreaction(Model, 'mRNA_deact + RNase <-> RNase_bound_d');
            RNase_binding_mRNAdeact_kineticlaw = addkineticlaw(RNase_binding_mRNAdeact,'MassAction');
            RNase_binding_mRNAdeact_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_F','TXTL_RNAdeg_R'};
            

                % Add degradation reaction with RNase_bound_d
            for stoich_coeff = 0:1:2
                RNase_bound_d_deg = addreaction(Model,'RNase_bound_d -> RNase');
                if ~isequal(stoich_coeff,0)
                    addproduct(RNase_bound_d_deg,mRNA_deact,stoich_coeff); 
                end
                RNase_bound_d_deg_kineticLaw = addkineticlaw(RNase_bound_d_deg,'MassAction'); 
                    % Modify kinetic parameter
                eval(sprintf('parameter_var_oi = ''TXTL_RNAdeg_kc_%s'';',num2str(stoich_coeff)))
                RNase_bound_d_deg.KineticLaw.ParameterVariableNames = {parameter_var_oi}; 
            end

            % Clean up - remove duplicate reactions, remove unused species
            rxn_idx = 1;
            while rxn_idx < length(Model.Reactions)
                rxn_oi = Model.Reactions(rxn_idx); 
                rxn_name = rxn_oi.Reaction; 
                all_model_rxn_name = {Model.Reactions.Reaction}; 
                rxn_name_idx_list = find(strcmp(rxn_name,all_model_rxn_name)); 
                if length(rxn_name_idx_list) > 1 
                    rxn_to_be_deleted = Model.Reactions(rxn_name_idx_list(2)); 
                    delete(rxn_to_be_deleted); 
                else
                    rxn_idx = rxn_idx + 1; 
                end
            end
            
               % Remove species that are not used in any reactions
            reaction_name_list = {Model.Reactions.Reaction};
            species_idx = 1; 
            while species_idx <= length(Model.Species)
                selected_species = Model.Species(species_idx); 
                exist_flag = false; 
                for rxn_idx = 1:length(reaction_name_list)
                    rxn_name = reaction_name_list{rxn_idx}; 
                    if contains(rxn_name,selected_species.Name)
                        exist_flag = true; 
                    end
                end
                if ~exist_flag
                    delete(selected_species)
                else
                    species_idx = species_idx + 1; 
                end
            end

            %% Track available binding sites for RNase 
        case 'TestBindingSites'

            % Collect all mRNA species in the model
            all_species_name_list = {Model.Species.Name};

            % Let's first add new species
                % Deactivated species should be agnostic to mRNA species
            RNase_bound_deactivated = addspecies(Model,'RNase_bound_d');
            binding_site_deactivated = addspecies(Model,'binding_site_d');

            % Identify all RNA species present 
            RNA_species_idx_list = contains(all_species_name_list,'RNA utr') & ~contains(all_species_name_list,':'); 
            RNA_species_name_list = all_species_name_list(RNA_species_idx_list); 

            for species_idx = 1:length(RNA_species_name_list)

                RNA_species_name = RNA_species_name_list{species_idx};
                RNA_species_name_for_model = strrep(RNA_species_name,' ','_');
                RNA_species_name_for_model = strrep(RNA_species_name_for_model,'-','');

                    % Add functional species name for reporter
                if contains(RNA_species_name,'sfGFP') || contains(RNA_species_name,'broc')
                    RNase_bound_functional_species_name = sprintf('RNase_bound_f_%s',RNA_species_name_for_model);
                    RNase_bound_functional = addspecies(Model,RNase_bound_functional_species_name);
                    binding_site_functional_species_name = sprintf('binding_site_f_%s',RNA_species_name_for_model);
                    binding_site_functional = addspecies(Model,binding_site_functional_species_name);
                end

                    % Add parameters 
                num_RNaseBindingSite_param_name = sprintf('num_RNaseBindingSite_%s',RNA_species_name_for_model);
                eval(sprintf('num_RNaseBindingSite_%s = addparameter(Model,num_RNaseBindingSite_param_name,5);',RNA_species_name_for_model)); 

                    % Assign initial binding site concentration 
                if contains(RNA_species_name,'sfGFP') || contains(RNA_species_name,'broc')
                    initial_site_functional_assignment_event_name = sprintf('binding_site_f_%s = [%s]',RNA_species_name_for_model,RNA_species_name);
                    initial_site_functional_assignment_event = addevent(Model,'time > 0',initial_site_functional_assignment_event_name); 

                    % Add binding reactions with sites 
                    RNase_binding_sitef_rxn_string = sprintf('binding_site_f_%s + RNase <-> RNase_bound_f_%s',RNA_species_name_for_model,RNA_species_name_for_model); 
                    RNase_binding_sitef = addreaction(Model,RNase_binding_sitef_rxn_string);
                    RNase_binding_sitef_kineticlaw = addkineticlaw(RNase_binding_sitef,'MassAction');
                    set(RNase_binding_sitef_kineticlaw,'ParameterVariableNames',{'TXTL_RNAdeg_F','TXTL_RNAdeg_R'}); 
                    % RNase_binding_sitef_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_F','TXTL_RNAdeg_R'};

                    % Add degradation reactions for the newly added
                    % RNase_bound_f species 
                    RNase_bound_f_deg_rxn_string = sprintf('RNase_bound_f_%s -> null',RNA_species_name_for_model); 
                    RNase_bound_f_deg_rxn = addreaction(Model,RNase_bound_f_deg_rxn_string); 
                    RNase_bound_f_deg_rxn_kineticLaw = addkineticlaw(RNase_bound_f_deg_rxn,'MassAction');
                    set(RNase_bound_f_deg_rxn_kineticLaw,'ParameterVariableNames',{'TXTL_RNAdeg_kc'});

                    % Add new reporter RNA -> null reactions with the same
                    % reaction rate 
                    deg_rxn_string = sprintf('%s -> null',RNA_species_name);
                    %         % find the binding rxn 
                    % all_model_rxn = {Model.Reactions.Reaction};
                    % num_occurrence_cell = strfind(all_model_rxn,reactant_RNA_name_for_model);
                    % num_occurrence = cellfun(@length,num_occurrence_cell,'UniformOutput',false);
                    % rxn_selected_for_mRNAdeg_rate = Model.Reactions(find(cell2mat(num_occurrence) > 1)).ReactionRate;
                    mRNA_degradation_rxn = addreaction(Model,deg_rxn_string,'ReactionRate',RNase_bound_f_deg_rxn.ReactionRate); 
                    mRNA_degradation_rxn_kineticLaw = addkineticlaw(mRNA_degradation_rxn,'Unknown'); 

                else % construct a long assignment string for all deactivated binding sites
                    if ~exist('initial_site_deactivated_assignment_string','var')
                        initial_site_deactivated_assignment_string = sprintf('binding_site_d = num_RNaseBindingSite_%s * [%s]',RNA_species_name_for_model,RNA_species_name);
                    else
                        initial_site_deactivated_assignment_string = strcat(initial_site_deactivated_assignment_string,sprintf(' + num_RNaseBindingSite_%s * [%s]'...
                            ,RNA_species_name_for_model,RNA_species_name));
                    end
                end
            end

            initial_site_deactivated_assignment_event = addevent(Model,'time > 0',initial_site_deactivated_assignment_string); 

            RNase_binding_sited = addreaction(Model, 'binding_site_d + RNase <-> RNase_bound_d');
            RNase_binding_sited_kineticlaw = addkineticlaw(RNase_binding_sited,'MassAction');
            RNase_binding_sited_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_F','TXTL_RNAdeg_R'};

            % Add degradation reaction (deactivated RNase_bound)
            RNase_deactivated_degradation = addreaction(Model, 'RNase_bound_d -> RNase');
            RNase_deactivated_degradation_kineticlaw = addkineticlaw(RNase_deactivated_degradation,'MassAction');
            RNase_deactivated_degradation_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_kc'};

            % Add degradation parameters to the model
            TXTL_RNAdeg_F = addparameter(Model,'TXTL_RNAdeg_F',1); 
            TXTL_RNAdeg_R = addparameter(Model,'TXTL_RNAdeg_R',1); 
            TXTL_RNAdeg_kc = addparameter(Model,'TXTL_RNAdeg_kc',1); 

            % Modify reactant 
            % Delete the original binding & degradation reaction
            rxn_to_be_deleted = {}; 
            for rxn_idx = 1:length(Model.Reactions)

                % Get both sides on the chemical equation 
                reaction_selected = Model.Reactions(rxn_idx); 
                chemical_reaction = reaction_selected.Reaction; 
                divider_idx = strfind(chemical_reaction,'->'); 
                rxn_lhs = chemical_reaction(1:divider_idx); 
                rxn_rhs = chemical_reaction(divider_idx + 1:end);
                rxn_rate = get(reaction_selected,'ReactionRate'); 
                rxn_kineticParams = reaction_selected.KineticLaw.Parameters; 
                % rxn_kineticParamNames = reaction_selected.KineticLaw.ParameterVariableNames;

                if contains(rxn_lhs,' RNase ') && contains(rxn_rhs,':RNase')
                    rxn_to_be_deleted{end + 1} = reaction_selected; 
                end

                if contains(rxn_lhs,':RNase') && contains(rxn_rhs,' RNase ')
                    rxn_to_be_deleted{end + 1} = reaction_selected; 
                end

                % Add binding site functional and binding site deactivated to any reaction that produces
                % mRNA 
                    % if RNA species list in rhs
                for RNA_species_name_idx = 1:length(RNA_species_name_list)

                    RNA_species_name = RNA_species_name_list{RNA_species_name_idx}; 
                    if contains(rxn_rhs,RNA_species_name) && contains(rxn_rhs,'term_') && ~contains(rxn_rhs,'RNase')

                        RNA_species_name = strrep(RNA_species_name,' ','_');
                        RNA_species_name = strrep(RNA_species_name,'-','');

                        if contains(RNA_species_name,'sfGFP') || contains(RNA_species_name,'broc')
                            product_name = sprintf('binding_site_f_%s',RNA_species_name); 
                            addproduct(reaction_selected,product_name,1)
                        end

                        % generation = sum(num_RNaseBindingSite * rate_of_mRNA production)
                        binding_site_deactivated_generation_rxn_rate_string = sprintf('num_RNaseBindingSite_%s * (%s)',RNA_species_name,rxn_rate);
                        % Add binding site generation reaction 
                        eval(sprintf('binding_site_deactivated_generation_%s%s = addreaction(Model, ''null -> binding_site_d'');',...
                            reaction_selected.Name,RNA_species_name));
                        eval(sprintf('kineticLaw_%s%s = addkineticlaw(binding_site_deactivated_generation_%s%s, "Unknown");',...
                            reaction_selected.Name,RNA_species_name,reaction_selected.Name,RNA_species_name)); 
                        eval(sprintf('set(binding_site_deactivated_generation_%s%s,''ReactionRate'',binding_site_deactivated_generation_rxn_rate_string);',...
                            reaction_selected.Name,RNA_species_name)); 
                            % Add kinetic parameters to the reaction 
                        for rxn_param_idx = 1:length(rxn_kineticParams)
                            rxn_kineticParam = rxn_kineticParams(rxn_param_idx); 
                            eval(sprintf('addparameter(kineticLaw_%s%s,rxn_kineticParam.Name,rxn_kineticParam.Value);',reaction_selected.Name,RNA_species_name));
                        end
                        % rxn_kineticParamNames = {rxn_kineticParams.Name}; 
                        % eval(sprintf('set(kineticLaw_%s%s, ''ParameterVariableNames'', rxn_kineticParamNames);',reaction_selected.Name,RNA_species_name)); 
                        % 
                    end

                end
                

            end

            

            for delete_rxn_idx = 1:length(rxn_to_be_deleted)
                delete_rxn_selected = rxn_to_be_deleted{delete_rxn_idx};
                delete(delete_rxn_selected)
            end
               

            % Clean up - remove duplicate reactions, remove unused species
            rxn_idx = 1;
            while rxn_idx < length(Model.Reactions)
                rxn_oi = Model.Reactions(rxn_idx); 
                rxn_name = rxn_oi.Reaction; 
                all_model_rxn_name = {Model.Reactions.Reaction}; 
                rxn_name_idx_list = find(strcmp(rxn_name,all_model_rxn_name)); 
                if length(rxn_name_idx_list) > 1 && ~contains(rxn_name,'null -> binding_site_d')
                    rxn_to_be_deleted = Model.Reactions(rxn_name_idx_list(2)); 
                    delete(rxn_to_be_deleted); 
                else
                    rxn_idx = rxn_idx + 1; 
                end
            end
            
               % Remove species that are not used in any reactions
            reaction_name_list = {Model.Reactions.Reaction};
            species_idx = 1; 
            while species_idx <= length(Model.Species)
                selected_species = Model.Species(species_idx); 
                exist_flag = false; 
                for rxn_idx = 1:length(reaction_name_list)
                    rxn_name = reaction_name_list{rxn_idx}; 
                    if contains(rxn_name,selected_species.Name)
                        exist_flag = true; 
                    end
                end
                if ~exist_flag
                    delete(selected_species)
                else
                    species_idx = species_idx + 1; 
                end
            end


    end
            % For newly added parameters/reactions, add parameters in Reaction.KineticLaw to model parameters and delete
        % them from kinetic law object 
    Mobj_rxns = get(Model,'Reactions');
    for rxn_idx = 1:length(Mobj_rxns)
        rxn = Mobj_rxns(rxn_idx);
        try
            rxn_params = rxn.KineticLaw.Parameters;
        catch
            rxn_params = []; 
        end
        for rxn_param_idx = 1:length(rxn_params)
            rxn_param = rxn_params(rxn_param_idx);
            rxn_param_name = rxn_param.Name; 
    
            if ~strcmp(rxn_param_name,rxn_param.Name)
                % Change parameter variable names
                parameter_name_idx = strcmp(rxn.KineticLaw.ParameterVariableNames,rxn_param.Name); 
                rxn.KineticLaw.ParameterVariableNames{parameter_name_idx} = rxn_param_name;
            end
            % Is this parameter already in model parameters?
            global_param = sbioselect(Model.Parameters,'Name',rxn_param_name);
            if ~isempty(global_param)
                % if already in model parameters, select and delete
                fprintf('\n %s already exists in model object,deleting',rxn_param_name)
                delete(rxn_param); 
            else
                % If doesn't exist in model parameters, move to model
                % parameters 
                    % Make sure parameter name matches 
                if ~strcmp(rxn_param_name,rxn_param.Name)
                    %If doesn't match, set parameter name
                    set(rxn_param,'Name',rxn_param_name); 
                end
                test_obj = move(rxn_param,Model);
            end
        end
        
    end







end