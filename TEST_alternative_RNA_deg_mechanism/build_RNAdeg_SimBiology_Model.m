function Model = build_RNAdeg_SimBiology_Model(kinetics)

    % Build (or Load) Simbiology model for RNA degradation 

    switch kinetics
        %% First-order kinetics RNA degradation 
        case 'firstOrder'
                % Create model object 
            Model = sbiomodel('m1');

                % Add species
            RNA = addspecies(Model,'RNA utrbroc--no_protein'); 
            AGMP = addspecies(Model,'AGMP'); 
            CUMP = addspecies(Model,'CUMP'); 

                % Add reaction and kinetic law 
            RNA_deg_1storder = addreaction(Model,'[RNA utrbroc--no_protein] -> 38 AGMP + 38 CUMP');
            kineticLawObj = addkineticlaw(RNA_deg_1storder,'MassAction');
            kineticParam = addparameter(kineticLawObj,'k_deg'); 
            kineticLawObj.ParameterVariableNames = 'k_deg';

        %% Michaelis-Menten kinetics RNA degradation 
        case 'MichaelisMenten'
                % Create model object 
            Model = sbiomodel('m1');

                % Add species
            RNA = addspecies(Model,'RNA utrbroc--no_protein'); 
            AGMP = addspecies(Model,'AGMP'); 
            CUMP = addspecies(Model,'CUMP'); 

                % Add reaction and kinetic law 
            RNA_deg_MM = addreaction(Model,'[RNA utrbroc--no_protein] -> 38 AGMP + 38 CUMP');
            kineticLawObj = addkineticlaw(RNA_deg_MM,'Henri-Michaelis-Menten');
            Vm = addparameter(kineticLawObj,'Vm'); 
            Km = addparameter(kineticLawObj,'Km'); 
            kineticLawObj.ParameterVariableNames = {'Vm','Km'};
            kineticLawObj.SpeciesVariableNames = 'RNA utrbroc--no_protein';

        %% Mass-action kinetics for RNA degradation 
        case 'MassAction'

            Model_struct = sbioloadproject('simbiology_models/RNA_deg.sbproj');
            Model = Model_struct.m1; 

        %% Accounting for deactivated mRNA species in degradation 
        case 'MixedDeg'

            Model_struct = sbioloadproject('simbiology_models/RNA_deg.sbproj');
            Model = Model_struct.m1; 

            % Let's first add new species
            RNase_bound = addspecies(Model,'RNase_bound');
            mRNA_deact = addspecies(Model,'mRNA_deact');
            
            % (1) For each RNase binding reaction, change the product into RNase_bound
            % (2) Duplicate the RNase degradation reaction and add deactivated mRNA as a product 
            reactions_to_be_added = {}; 
            for rxn_idx = 1:length(Model.Reactions)
            
                reaction_selected = Model.Reactions(rxn_idx); 
                chemical_reaction = reaction_selected.Reaction; 
                divider_idx = strfind(chemical_reaction,'->'); 
                rxn_lhs = chemical_reaction(1:divider_idx); 
                rxn_rhs = chemical_reaction(divider_idx + 1:end); 
                % find RNase binding reactions and replace product 
                if contains(rxn_lhs,' RNase ') && contains(rxn_rhs,':RNase')  
            
                    % Get product name and remove product 
                    product_name = get(reaction_selected.product); 
                    rmproduct(reaction_selected, product_name.Name); 
            
                    % Replace with RNase_bound
                    addproduct(reaction_selected,RNase_bound,1); 
                end
            
                % Find RNA degradation reaction and (1) replace reactant (2) get ready to duplicate those 
                % (3) Remove AGMP/CUMP 
                if contains(rxn_lhs,':RNase') && contains(rxn_rhs,' RNase ')
            
                    reactions_to_be_added{end + 1} = reaction_selected;
            
                    % Get reactant name and remove reactant 
                    reactant_name = get(reaction_selected.Reactants);
                    rmreactant(reaction_selected,reactant_name.Name); 
            
                    % Replace with RNase_bound
                    addreactant(reaction_selected,RNase_bound,1); 
            
                    % Remove AGMP/CUMP 
                    rmproduct(reaction_selected,'AGMP');
                    rmproduct(reaction_selected,'CUMP'); 
            
                end
            
            end
            
            % Add a new reaction to account for RNase binding to deactivated mRNA 
                % Add paramet
            TXTL_RNAdeg_R_proc = addparameter(Model,'TXTL_RNAdeg_kc_proc',1); 
            
              % Add new reactions
            RNase_binding_mRNAdeact = addreaction(Model, 'mRNA_deact + RNase <-> RNase_bound');
            RNase_binding_mRNAdeact_kineticlaw = addkineticlaw(RNase_binding_mRNAdeact,'MassAction');
            RNase_binding_mRNAdeact_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_F','TXTL_RNAdeg_R'};
            
            % Duplciate the reactions and add deactivated mRNA as a product 
            for dup_rxn_idx = 1:length(reactions_to_be_added)
                rxn_oi = reactions_to_be_added{dup_rxn_idx};
                dup_rxn_oi = copyobj(rxn_oi,Model); 
            
                % Add product to the duplicated reaction 
                addproduct(dup_rxn_oi,mRNA_deact,1); 
            
            end
            
            % Use different kinetic parameters for processive and non-processive RNA
            % degradation 
                % These should act on the same set of reactions 
            for proc_deg_rxn_idx = 1:length(reactions_to_be_added)
                rxn_oi = reactions_to_be_added{proc_deg_rxn_idx};
            
                % Change the associated parameter variable names 
                rxn_oi.KineticLaw.ParameterVariableNames = {'TXTL_RNAdeg_kc_proc'}; 
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

        %% Tracking functional & deactivated mRNA-RNase complexes separately 
        case 'TestPrototype'

            Model_struct = sbioloadproject('simbiology_models/RNA_deg.sbproj');
            Model = Model_struct.m1; 

            % Let's first add new species
            RNase_bound_functional = addspecies(Model,'RNase_bound_f');
            RNase_bound_deactivated = addspecies(Model,'RNase_bound_d');
            mRNA_deact = addspecies(Model,'mRNA_deact');

            % (1) For each RNase binding reaction, change the product into
            % RNase_bound_f
            % (2) Record these reactions to add 3 reactions later to each
            binding_rxns_to_be_added = {}; 
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

                    % Add to the tracking pool
                    binding_rxns_to_be_added{end + 1} = reaction_selected;
            
                    % Get product name and remove product 
                    product_name = get(reaction_selected.product); 
                    rmproduct(reaction_selected, product_name.Name); 
            
                    % Replace with RNase_bound
                    addproduct(reaction_selected,RNase_bound_functional,1); 
                end
            
                % Find RNA degradation reaction and (1) replace reactant
                % (2) add to tracking 
                % (3) Remove AGMP/CUMP (for now) 
                if contains(rxn_lhs,':RNase') && contains(rxn_rhs,' RNase ')
            
                    deg_rxns_to_be_added{end + 1} = reaction_selected;
            
                    % Get reactant name and remove reactant 
                    reactant_name = get(reaction_selected.Reactants);
                    rmreactant(reaction_selected,reactant_name.Name); 
            
                    % Replace with RNase_bound
                    addreactant(reaction_selected,RNase_bound_functional,1); 

                    % Add 2 deactivated mRNA as product 
                    addproduct(reaction_selected,mRNA_deact,2)
            
                    % Remove AGMP/CUMP 
                    rmproduct(reaction_selected,'AGMP');
                    rmproduct(reaction_selected,'CUMP'); 
            
                end
            
            end

                % Add parameters to represent degradation kinetic constant
                % for deactivated mRNAs; these parameters are technically
                % time-dependent but the average will be taken 
            TXTL_RNAdeg_kc_2 = addparameter(Model,'TXTL_RNAdeg_kc_2',1); 
            TXTL_RNAdeg_kc_1 = addparameter(Model,'TXTL_RNAdeg_kc_1',1); 
            TXTL_RNAdeg_kc_0 = addparameter(Model,'TXTL_RNAdeg_kc_0',1); 
            
              % Add binding reaction with mRNA_deactivated 
            RNase_binding_mRNAdeact = addreaction(Model, 'mRNA_deact + RNase <-> RNase_bound_d');
            RNase_binding_mRNAdeact_kineticlaw = addkineticlaw(RNase_binding_mRNAdeact,'MassAction');
            RNase_binding_mRNAdeact_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_F','TXTL_RNAdeg_R'};
            
            % Duplciate the degradation reactions and add deactivated mRNA as a product 
            for dup_rxn_idx = 1:length(deg_rxns_to_be_added)
                % Add product to the duplicated reaction 
                for stoich_coeff = 0:1:2
                    rxn_oi = deg_rxns_to_be_added{dup_rxn_idx};
                    dup_rxn_oi = copyobj(rxn_oi,Model); 
                    
                    % Replace reactant with RNase_bound_d
                        % Get reactant name and remove reactant 
                    reactant_name = get(dup_rxn_oi.Reactants);
                    rmreactant(dup_rxn_oi,reactant_name.Name); 
            
                        % Replace with RNase_bound
                    addreactant(dup_rxn_oi,RNase_bound_deactivated,1); 
            
                    % Remove existing deactivated mRNA & Add deactivated mRNA to the product side 
                    rmproduct(dup_rxn_oi,'mRNA_deact'); 
                    if ~isequal(stoich_coeff,0)
                        addproduct(dup_rxn_oi,mRNA_deact,stoich_coeff); 

                    end

                    % Modify kinetic parameter
                    eval(sprintf('parameter_var_oi = ''TXTL_RNAdeg_kc_%s'';',num2str(stoich_coeff)))
                    dup_rxn_oi.KineticLaw.ParameterVariableNames = {parameter_var_oi}; 
                end

            
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

            Model_struct = sbioloadproject('simbiology_models/RNA_deg.sbproj');
            Model = Model_struct.m1; 

            % Add an initial assignment rule to account for initial site
            % concentration 
                % Collect all mRNA species in the model
            all_species_name_list = {Model.Species.Name};

            % Let's first add new species
                % Deactivated species should be agnostic to mRNA species
            RNase_bound_deactivated = addspecies(Model,'RNase_bound_d');
            binding_site_deactivated = addspecies(Model,'binding_site_d');

            % Delete the original binding reaction
            rxn_to_be_deleted = {}; 
            for rxn_idx_2 = 1:length(Model.Reactions)
            
                % Get both sides on the chemical equation 
                reaction_selected = Model.Reactions(rxn_idx_2); 
                chemical_reaction = reaction_selected.Reaction; 
                divider_idx = strfind(chemical_reaction,'->'); 
                rxn_lhs = chemical_reaction(1:divider_idx); 
                rxn_rhs = chemical_reaction(divider_idx + 1:end);
                if contains(rxn_lhs,' RNase ') && contains(rxn_rhs,':RNase')
                    rxn_to_be_deleted{end + 1} = reaction_selected; 
                end

            end
            for rxn_delete_idx = 1:length(rxn_to_be_deleted)
                rxn_selected = rxn_to_be_deleted{rxn_delete_idx}; 
                delete(rxn_selected)
            end

            % Identify all RNA species present 
            RNA_species_idx_list = contains(all_species_name_list,'RNA utr') & ~contains(all_species_name_list,':'); 
            RNA_species_name_list = all_species_name_list(RNA_species_idx_list); 

            for species_idx = 1:length(RNA_species_name_list)
                RNA_species_name = RNA_species_name_list{species_idx};
                RNA_species_name_for_model = strrep(RNA_species_name,' ','_');
                RNA_species_name_for_model = strrep(RNA_species_name_for_model,'-','');

                    % Add functional species name
                RNase_bound_functional_species_name = sprintf('RNase_bound_f_%s',RNA_species_name_for_model);
                RNase_bound_functional = addspecies(Model,RNase_bound_functional_species_name);
                binding_site_functional_species_name = sprintf('binding_site_f_%s',RNA_species_name_for_model);
                binding_site_functional = addspecies(Model,binding_site_functional_species_name);

                    % Add parameters 
                num_RNaseBindingSite_param_name = sprintf('num_RNaseBindingSite_%s',RNA_species_name_for_model);
                eval(sprintf('num_RNaseBindingSite_%s = addparameter(Model,num_RNaseBindingSite_param_name,5);',RNA_species_name_for_model)); 


                initial_site_functional_assignment_event_name = sprintf('binding_site_f_%s = [%s]',RNA_species_name_for_model,RNA_species_name);
                initial_site_functional_assignment_event = addevent(Model,'time > 0',initial_site_functional_assignment_event_name); 
        
                initial_site_deactivated_assignment_event_name = sprintf('binding_site_d = num_RNaseBindingSite_%s * [%s]',...
                        RNA_species_name_for_model,RNA_species_name);
                initial_site_deactivated_assignment_event = addevent(Model,'time > 0',initial_site_deactivated_assignment_event_name); 

                % Add binding reactions with sites 
                RNase_binding_sitef_rxn_string = sprintf('binding_site_f_%s + RNase <-> RNase_bound_f_%s',RNA_species_name_for_model,RNA_species_name_for_model); 
                RNase_binding_sitef = addreaction(Model,RNase_binding_sitef_rxn_string);
                RNase_binding_sitef_kineticlaw = addkineticlaw(RNase_binding_sitef,'MassAction');
                RNase_binding_sitef_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_F','TXTL_RNAdeg_R'};
            end

            RNase_binding_sited = addreaction(Model, 'binding_site_d + RNase <-> RNase_bound_d');
            RNase_binding_sited_kineticlaw = addkineticlaw(RNase_binding_sited,'MassAction');
            RNase_binding_sited_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_F','TXTL_RNAdeg_R'};

            % Add degradation reaction (deactivated RNase_bound)
            RNase_deactivated_degradation = addreaction(Model, 'RNase_bound_d -> RNase');
            RNase_deactivated_degradation_kineticlaw = addkineticlaw(RNase_deactivated_degradation,'MassAction');
            RNase_deactivated_degradation_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_kc'};

            % Modify reactant 
            for rxn_idx = 1:length(Model.Reactions)
            
                % Get both sides on the chemical equation 
                reaction_selected = Model.Reactions(rxn_idx); 
                chemical_reaction = reaction_selected.Reaction; 
                divider_idx = strfind(chemical_reaction,'->'); 
                rxn_lhs = chemical_reaction(1:divider_idx); 
                rxn_rhs = chemical_reaction(divider_idx + 1:end);
                rxn_rate = get(reaction_selected,'ReactionRate'); 

                if contains(rxn_lhs,':RNase') && contains(rxn_rhs,' RNase ')
            
                    % Get reactant name and remove reactant 
                    reactant_name = get(reaction_selected.Reactants);
                    reactant_RNA_name = reactant_name.Name(1:strfind(reactant_name.Name,':')-1);
                    rmreactant(reaction_selected,reactant_name.Name); 
            
                    % Replace with RNase_bound_functional
                    reactant_RNA_name_for_model = strrep(reactant_RNA_name,' ','_');
                    reactant_RNA_name_for_model = strrep(reactant_RNA_name_for_model,'-','');
                    eval(sprintf('addreactant(reaction_selected,"RNase_bound_f_%s",1);',reactant_RNA_name_for_model)); 
            
                    % Remove AGMP/CUMP 
                    rmproduct(reaction_selected,'AGMP');
                    rmproduct(reaction_selected,'CUMP'); 

                    % Also add a parallel degradation reaction for mRNA
                    % with the same reaction kinetics 

                        % Add reaction with assigned rate 
                    deg_rxn_string = sprintf('%s -> null',reactant_RNA_name);
                    %         % find the binding rxn 
                    % all_model_rxn = {Model.Reactions.Reaction};
                    % num_occurrence_cell = strfind(all_model_rxn,reactant_RNA_name_for_model);
                    % num_occurrence = cellfun(@length,num_occurrence_cell,'UniformOutput',false);
                    % rxn_selected_for_mRNAdeg_rate = Model.Reactions(find(cell2mat(num_occurrence) > 1)).ReactionRate;
                    mRNA_degradation_rxn = addreaction(Model,deg_rxn_string,'ReactionRate',reaction_selected.ReactionRate); 
            
                end

                % Add binding site functional and binding site deactivated to any reaction that produces
                % mRNA 
                    % if RNA species list in rhs
                for RNA_species_name_idx = 1:length(RNA_species_name_list)
                    RNA_species_name = RNA_species_name_list{RNA_species_name_idx}; 
                    if contains(rxn_rhs,RNA_species_name) && ~contains(rxn_rhs,':')
                        RNA_species_name = strrep(RNA_species_name,' ','_');
                        RNA_species_name = strrep(RNA_species_name,'-','');
                        product_name = sprintf('binding_site_f_%s',RNA_species_name); 
                        addproduct(reaction_selected,product_name,1)
                        % generation = num_RNaseBindingSite * rate_of_mRNA production 
                        binding_site_deactivated_generation_rxn_rate_string = sprintf('num_RNaseBindingSite_%s * (%s)',RNA_species_name,rxn_rate);
                        binding_site_deactivated_generation = addreaction(Model, 'null -> binding_site_d',...
                            'ReactionRate',binding_site_deactivated_generation_rxn_rate_string);
                        
                    end

                end

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