function problemObject = setProblemObject_v3(group_description,group_number,output_path)
    % Construct a problemObject containing (1) the master model for plasmid
    % crosstalk modeling (2) relevant data required for parameter
    % estimation, as defined by group_description (3) information on
    % parameter to be estimated, including initial values, LB & UB 

    % This is designed to fit both transcriptional-level and protein-level crosstalk.

    % var3 assumes that translation generates the toxin and it acts on
    % ribosomes

    % INPUT:
        % GROUP DESCRIPTION: Option to select experimental data to be fitted by bundling all plasmid
        % concentrations; a cell array of text descriptions 
        % GROUP NUMBER: Option to select experimental data by corresponding group
        % numbers and TX/protein; a dictionary with 'TX' or 'PE' as keys.
        % Values correpond the data group number in each case (1-84 for TX,
        % 1-112 for PE) 
        % OUTPUT PATH: path to parameter initial values, LB & UB

        % Note: Either group description or group number should be used.
        % While one is used, the other one should be an empty vector 

    % v3 update: (1) check model for duplicate reactions and remove those
    % (2) Remove weights, error model, and short time course options (3)
    % Add group number option to include only selected data 

    %% Set correct path to avoid confusion of config files
    currentpath = pwd;
    addpath(genpath(sprintf('%s/plasmid_crosstalk_config_files',currentpath))); 
    addpath(genpath(sprintf('%s/txtlsim_buildacell',currentpath))); 
    addpath(genpath(sprintf('%s',currentpath)));
    rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath));
    problemObject = fitproblem(); % Initialize a fitProblem object
    
    %% Load data 

        % load both transcription and protein expression data 
    tx_data_file = load('data_structures/simbio_data_table_updated_FP.mat'); 
    PE_data_file = load('data_structures/simbio_data_table_PE_updated_FP.mat');
    tx_data_table = tx_data_file.all_data_table; 
    PE_data_table = PE_data_file.all_data_table; 
    tx_data_description = tx_data_file.variable_name_stem_list; 
    PE_data_description = PE_data_file.variable_name_stem_list; 

        % Select relevant data 
    if ~isempty(group_description) && ~isempty(group_number)
        error('Both data group description and group number are provided. Only a single option is allowed')
    elseif isempty(group_description) && isempty(group_number)
        error('No data group description or group number given')
    elseif ~isempty(group_description)
        [~,tx_selected_table] = get_relevant_data(tx_data_table,tx_data_description,group_description); 
        [~,PE_selected_table] = get_relevant_data(PE_data_table,PE_data_description,group_description);
    elseif ~isempty(group_number)
        % TX
        TX_group_number = group_number({'TX'}); 
        TX_group_number = TX_group_number{1}; 
        tx_selected_idx = []; 
        for group_num = TX_group_number
            tx_selected_idx = [tx_selected_idx;find(tx_data_table.Group==group_num)];
        end
        tx_selected_table = tx_data_table(tx_selected_idx,:); 
    
                % PE 
        PE_group_number = group_number({'PE'}); 
        PE_group_number = PE_group_number{1}; 
        PE_selected_idx = []; 
        for group_num = PE_group_number
            PE_selected_idx = [PE_selected_idx;find(PE_data_table.Group==group_num)];
        end
        PE_selected_table = PE_data_table(PE_selected_idx,:); 
    end

        % Convert fluorescence to concentration and add concentration to
        % data table 
        %%%%%% mRNA %%%%%%
    if ~isempty(tx_selected_table)
        split_point_idx_list = find(tx_selected_table.gain==tx_selected_table.gain(1)); 
        split_point_idx = split_point_idx_list(end); 
        if isequal(length(split_point_idx_list),height(tx_selected_table))% if there's only one gain 
            mRNA_concentration = [fluo_mRNA_conversion(tx_selected_table.fluorescence,[],tx_selected_table.gain(1))];
        else
                % Applies to when there are 2 gains 
            mRNA_concentration = [fluo_mRNA_conversion(tx_selected_table.fluorescence(1:split_point_idx),[],tx_selected_table.gain(1));...
                fluo_mRNA_conversion(tx_selected_table.fluorescence(split_point_idx+1:end),[],tx_selected_table.gain(split_point_idx + 1))];
        end
        tx_selected_table = addvars(tx_selected_table,mRNA_concentration,'NewVariableNames','mRNA_concentration'); 
            % Replace negative values with nan 
        tx_selected_table = remove_neg_val(tx_selected_table);
            % Find new split after negative value removal
        new_split_point_idx_list = find(tx_selected_table.gain==tx_selected_table.gain(1)); 
        new_split_point_idx = new_split_point_idx_list(end); 
        if isequal(length(new_split_point_idx_list),height(tx_selected_table))
            mRNA_CI_lb = [fluo_mRNA_conversion(tx_selected_table.CI_lb,[],tx_selected_table.gain(1))];
            mRNA_CI_ub = [fluo_mRNA_conversion(tx_selected_table.CI_ub,[],tx_selected_table.gain(1))];
        else
            mRNA_CI_lb = [fluo_mRNA_conversion(tx_selected_table.CI_lb(1:new_split_point_idx),[],tx_selected_table.gain(1));...
                fluo_mRNA_conversion(tx_selected_table.CI_lb(new_split_point_idx+1:end),[],tx_selected_table.gain(new_split_point_idx + 1))];
            mRNA_CI_ub = [fluo_mRNA_conversion(tx_selected_table.CI_ub(1:new_split_point_idx),[],tx_selected_table.gain(1));...
                fluo_mRNA_conversion(tx_selected_table.CI_ub(new_split_point_idx+1:end),[],tx_selected_table.gain(new_split_point_idx + 1))];
        end
        tx_selected_table = addvars(tx_selected_table,mRNA_CI_lb,'NewVariableNames','mRNA_concentration_CI_lb'); 
        tx_selected_table = addvars(tx_selected_table,mRNA_CI_ub,'NewVariableNames','mRNA_concentration_CI_ub'); 
    end
    %%%%%% Protein %%%%%%
    if ~isempty(PE_selected_table)
        PE_split_point_idx_list = find(PE_selected_table.gain==PE_selected_table.gain(1)); 
        PE_split_point_idx = PE_split_point_idx_list(end); 
        if isequal(length(PE_split_point_idx_list),height(PE_selected_table))% if there's only one gain 
            GFP_concentration = [fluo_GFP_conversion(PE_selected_table.fluorescence,[],PE_selected_table.gain(1))];
        else
                % Applies to when there are 2 gains 
            GFP_concentration = [fluo_GFP_conversion(PE_selected_table.fluorescence(1:PE_split_point_idx),[],PE_selected_table.gain(1));...
                fluo_GFP_conversion(PE_selected_table.fluorescence(PE_split_point_idx+1:end),[],PE_selected_table.gain(PE_split_point_idx + 1))];
        end
        PE_selected_table = addvars(PE_selected_table,GFP_concentration,'NewVariableNames','GFP_concentration'); 
            % Replace negative values with nan 
        PE_selected_table = remove_neg_val(PE_selected_table);
            % Find new split after negative value removal
        new_split_point_idx_list = find(PE_selected_table.gain==PE_selected_table.gain(1)); 
        new_split_point_idx = new_split_point_idx_list(end); 
        if isequal(length(new_split_point_idx_list),height(PE_selected_table))
            GFP_CI_lb = [fluo_GFP_conversion(PE_selected_table.CI_lb,[],PE_selected_table.gain(1))];
            GFP_CI_ub = [fluo_GFP_conversion(PE_selected_table.CI_ub,[],PE_selected_table.gain(1))];
        else
            GFP_CI_lb = [fluo_GFP_conversion(PE_selected_table.CI_lb(1:new_split_point_idx),[],PE_selected_table.gain(1));...
                fluo_GFP_conversion(PE_selected_table.CI_lb(new_split_point_idx+1:end),[],PE_selected_table.gain(new_split_point_idx + 1))];
            GFP_CI_ub = [fluo_GFP_conversion(PE_selected_table.CI_ub(1:new_split_point_idx),[],PE_selected_table.gain(1));...
                fluo_GFP_conversion(PE_selected_table.CI_ub(new_split_point_idx+1:end),[],PE_selected_table.gain(new_split_point_idx + 1))];
        end
        PE_selected_table = addvars(PE_selected_table,GFP_CI_lb,'NewVariableNames','GFP_concentration_CI_lb'); 
        PE_selected_table = addvars(PE_selected_table,GFP_CI_ub,'NewVariableNames','GFP_concentration_CI_ub'); 
    end
       % Consolidate columns/variables in PE_selected_table and
       % tx_selected_table
    tx_null_var = zeros(height(tx_selected_table),1);
    PE_null_var = zeros(height(PE_selected_table),1); 
    
    tx_selected_table_vars = tx_selected_table.Properties.VariableNames; 
    PE_selected_table_vars = PE_selected_table.Properties.VariableNames;
    all_table_vars = [tx_selected_table_vars,PE_selected_table_vars];
    all_table_vars = unique(all_table_vars);

    for var_idx = 1:length(all_table_vars)
        var_name = all_table_vars{var_idx};
        if ~any(strcmp(var_name,tx_selected_table_vars))
            tx_selected_table  = addvars(tx_selected_table,tx_null_var,'NewVariableNames',var_name);
        end
        if ~any(strcmp(var_name,PE_selected_table_vars))
            PE_selected_table  = addvars(PE_selected_table,PE_null_var,'NewVariableNames',var_name);
        end
    end

    % Reindex group number before combining 
    num_group_tx_table = length(tx_data_description);
    PE_selected_table.Group = PE_selected_table.Group + num_group_tx_table; 
    selected_table = [tx_selected_table;PE_selected_table];

    % Add to problemObject
    problemObject.Data = groupedData(selected_table); % define a groupedData object 
    problemObject.Data.Properties.IndependentVariableName = 'Time'; 
    problemObject.Data.Properties.GroupVariableName = 'Group'; 

    %% Build model
    % Define components 
    tube1 = txtl_extract('plasmid_crosstalk');
    tube2 = txtl_buffer('plasmid_crosstalk');
    tube3 = txtl_newtube('gene_expression');

    % Add DNA template 
        % 3WJdB
    dna_broc_T7_strong = txtl_add_dna(tube3, 'pT7(23)', 'utrbroc(152)', 'no_protein',0, 'plasmid');	
    dna_broc_T7_weak = txtl_add_dna(tube3, 'pT773(23)', 'utrbroc(152)', 'no_protein',0, 'plasmid');	
    dna_broc_sigma70_strong = txtl_add_dna(tube3, 'pJ23119(35)', 'utrbroc(800)', 'no_protein',0, 'plasmid');
        % sfGFP
    dna_GFP_T7_strong = txtl_add_dna(tube3,'pT7(23)', 'utrGFP(57)', 'sfGFP(723)',0,'plasmid'); 
    dna_GFP_T7_weak = txtl_add_dna(tube3,'pT773(23)', 'utrGFP(57)', 'sfGFP(723)',0,'plasmid'); 
    dna_GFP_sigma70_strong = txtl_add_dna(tube3,'pJ23119(35)', 'utrGFP(57)', 'sfGFP(723)',0,'plasmid'); 
    dna_GFP_sigma70_weak = txtl_add_dna(tube3,'pJ23105(35)', 'utrGFP(57)', 'sfGFP(723)',0,'plasmid'); 
        % empty plasmids 
    dna_kan = txtl_add_dna(tube3,'pkanR(20)','utrkanR(33)','kanR(783)',0,'plasmid');
    dna_empty_T7 = txtl_add_dna(tube3,'pT7(23)','utrempty(51)','no_protein',0,'plasmid'); 
    dna_empty_sigma70 = txtl_add_dna(tube3,'pJ23119(35)','utrempty(51)','no_protein',0,'plasmid'); 

    Mobj = txtl_combine([tube1, tube2, tube3]);
    % txtl_runsim mode initializes reactions and parameters to model object 

    %%%%%set the reporter species initial concentration to be non-zero%%%%%%%
    if strcmp(problemObject.ErrorModel,'proportional')
        tx_reporter_species = sbioselect(Mobj,'Name','RNA utrbroc--no_protein');
        PE_reporter_species = sbioselect(Mobj,'Name','protein sfGFP*');
        set(tx_reporter_species,'InitialAmount',1e-08); 
        set(PE_reporter_species,'InitialAmount',1e-08); 
    end
    %%%%%%%%%%%%

    [initial_simData] = txtl_runsim(Mobj,14*60*60);

    % Check for duplicate reactions and remove them 
    rxn_idx = 1;
    while rxn_idx < length(Mobj.Reactions)
        rxn_oi = Mobj.Reactions(rxn_idx); 
        rxn_name = rxn_oi.Reaction; 
        all_model_rxn_name = {Mobj.Reactions.Reaction}; 
        rxn_name_idx_list = find(strcmp(rxn_name,all_model_rxn_name)); 
        if length(rxn_name_idx_list) > 1 
            rxn_to_be_deleted = Mobj.Reactions(rxn_name_idx_list(1)); 
            delete(rxn_to_be_deleted); 
        else
            rxn_idx = rxn_idx + 1; 
        end
    end

    %Add fields to problemObject
    problemObject.Model = Mobj; 

    %% Define the fitProblem struct 
    
        % % Change error model 
        % problemObject.ErrorModel = 'combined'; 

    % Link species in model to variable name of data 
        % Note that mRNA and GFP should be orthogonal (i.e. [RNA] = 0 or
        % [GFP] = 0) 
    if isempty(PE_selected_table)
        problemObject.ResponseMap = '[RNA utrbroc--no_protein] = mRNA_concentration';
    elseif isempty(tx_selected_table)
        problemObject.ResponseMap = '[protein sfGFP*] = GFP_concentration';
    else
        problemObject.ResponseMap = {'[RNA utrbroc--no_protein] = mRNA_concentration',...
            '[protein sfGFP*] = GFP_concentration'};
    end

    % Create doses from data
    problemObject.Doses = createDoses(problemObject.Data,{'T7_strong_reporter_plasmid','T7_weak_reporter_plasmid',...
        'sigma70_strong_reporter_plasmid','kan_plasmid','empty_T7_plasmid','empty_sigma70_plasmid',...
        'T7_strong_GFP_plasmid','T7_weak_GFP_plasmid','sigma70_strong_GFP_plasmid',...
        'sigma70_weak_GFP_plasmid'}); 
    for dose_idx = 1:size(problemObject.Doses,1)
        % This needs to be tailored to datasets used 
        problemObject.Doses(dose_idx,1).TargetName = 'DNA pT7--utrbroc--no_protein';
        problemObject.Doses(dose_idx,2).TargetName = 'DNA pT773--utrbroc--no_protein';
        problemObject.Doses(dose_idx,3).TargetName = 'DNA pJ23119--utrbroc--no_protein';
        problemObject.Doses(dose_idx,4).TargetName = 'DNA pkanR--utrkanR--kanR';
        problemObject.Doses(dose_idx,5).TargetName = 'DNA pT7--utrempty--no_protein'; 
        problemObject.Doses(dose_idx,6).TargetName = 'DNA pJ23119--utrempty--no_protein'; 
        problemObject.Doses(dose_idx,7).TargetName = 'DNA pT7--utrGFP--sfGFP'; 
        problemObject.Doses(dose_idx,8).TargetName = 'DNA pT773--utrGFP--sfGFP'; 
        problemObject.Doses(dose_idx,9).TargetName = 'DNA pJ23119--utrGFP--sfGFP'; 
        problemObject.Doses(dose_idx,10).TargetName = 'DNA pJ23105--utrGFP--sfGFP'; 
    end

    %% Estimated parameters info
    
    % Add parameters in Reaction.KineticLaw to model parameters and delete
    % them from kinetic law object 
    Mobj_rxns = get(Mobj,'Reactions');
    for rxn_idx = 1:length(Mobj_rxns)
        rxn = Mobj_rxns(rxn_idx);
        rxn_params = rxn.KineticLaw.Parameters;
        for rxn_param_idx = 1:length(rxn_params)
            rxn_param = rxn_params(rxn_param_idx);
            rxn_param_name = rxn_param.Name; 
            % If is RNA degradation parameters, modify name. 
            % Also modify ParameterVariableNames (Reaction rate
            % automatically updates) 
            if contains(rxn_param_name,'TXTL_RNAdeg_R') && contains(rxn.Reaction,'utrbroc')
                rxn_param_name = [rxn_param_name '_broc'];
            end
            if contains(rxn_param_name,'TXTL_RNAdeg_R') && (contains(rxn.Reaction,'sfGFP') || contains(rxn.Reaction,'empty'))
                rxn_param_name = [rxn_param_name '_sfGFP'];
            end
                % For now let's treat kanR and sfGFP as having different
                % degradation kinetic param
            if contains(rxn_param_name,'TXTL_RNAdeg_kc') && contains(rxn.Reaction,'utrkanR')
                rxn_param_name = [rxn_param_name '_kanR']; 
            end
            if contains(rxn_param_name,'TXTL_RNAdeg_kc') && contains(rxn.Reaction,'sfGFP')
                rxn_param_name = [rxn_param_name '_sfGFP']; 
            end
            if contains(rxn_param_name,'TXTL_RNAdeg_kc') && contains(rxn.Reaction,'broc')
                rxn_param_name = [rxn_param_name '_broc']; 
            end
            if ~strcmp(rxn_param_name,rxn_param.Name)
                % Change parameter variable names
                parameter_name_idx = strcmp(rxn.KineticLaw.ParameterVariableNames,rxn_param.Name); 
                rxn.KineticLaw.ParameterVariableNames{parameter_name_idx} = rxn_param_name;
            end
            % Is this parameter already in model parameters?
            global_param = sbioselect(Mobj.Parameters,'Name',rxn_param_name);
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
                test_obj = move(rxn_param,Mobj);
            end
        end
    end

    % Add polymerase, ribosome, and RNase initial concentrations as parameters
    % and add rules to set initial conditions for the species 
    species_init_unknown = {'RNAP','t7RNAP','RNase','Ribo'};
    for add_param_idx = 1:length(species_init_unknown)
        name_species_init_unknown = species_init_unknown{add_param_idx}; 
        % select the species 
        selected_species_obj = sbioselect(Mobj,'Type','Species','Name',name_species_init_unknown);
        % Add a parameter to model object with the name species.Name + _0 
        addparameter(Mobj,sprintf('%s_0',name_species_init_unknown),selected_species_obj.Value); 
        % Assign a rule to assign species initial value with parameter values 
        addrule(Mobj,sprintf('%s = %s_0',name_species_init_unknown,name_species_init_unknown)); 
    end

    % Taken transcription toxicity into account 
        % Add new parameters named 
    tx_capacity_param = addparameter(Mobj,'tx_capacity_param'); 
    k_toxin = addparameter(Mobj,'k_toxin'); 
    toxin_threshold = addparameter(Mobj,'toxin_threshold'); 
        % Add the rule to Mobj, both initial and repeated
        % assignments 
            % Put all RNA species into an array, if a reaction produces any
            % of those speicies, it should be added to the sum of
            % transcription rate 
           
            
    % model_species_name = get(Mobj.Species,'Name');
    % RNA_species_idx = contains(model_species_name,'RNA utr') & ~contains(model_species_name,'RNase') & ~contains(model_species_name,'Ribo'); 
    % RNA_species_name = model_species_name(RNA_species_idx); 
    % model_reaction_name = get(Mobj.Reactions,'Reaction'); 
    % tx_reaction_name_list = {};
    % for rxn_name_idx = 1:length(model_reaction_name)
    %     reaction_name_single = model_reaction_name{rxn_name_idx}; 
    %     rhs_start_idx = strfind(reaction_name_single,'->');
    %     truncated_test_reaction_name_single = reaction_name_single(rhs_start_idx:end); 
    %     for RNA_species_idx = 1:length(RNA_species_name)
    %         RNA_species_name_single = RNA_species_name{RNA_species_idx}; 
    %         if contains(truncated_test_reaction_name_single,strcat('[',RNA_species_name_single,']')) &&...
    %                 ~contains(truncated_test_reaction_name_single,'Ribo')
    %             tx_reaction_name_list{end+1} = reaction_name_single;
    %         end
    %     end
    % 
    % end

            % Nov 2023 edit: Instead of using transcription termination as the step to
            % produce toxin, add toxin to the elongation step (aka
            % consumption described in paper) 
    model_reaction_name_list = get(Mobj.Reactions,'Reaction');
    tl_elong_rxn_name_list = {};
    for rxn_name_idx = 1:length(model_reaction_name_list)
        reaction_name = model_reaction_name_list{rxn_name_idx}; 
        rhs_start_idx = strfind(reaction_name,'->');
        rxn_name_lhs = reaction_name(1:rhs_start_idx); 
        rxn_name_rhs = reaction_name(rhs_start_idx + 1:end); 
        if contains(rxn_name_lhs,'AA:AGTP:') && ~contains(rxn_name_rhs,'term') && ~contains(reaction_name,'RNase')
            tl_elong_rxn_name_list{end + 1} = reaction_name;
        end
    end

        % Add a null species as the product of these reactions 
    toxin = addspecies(Mobj,'toxin'); 
    set(toxin,'InitialAmount',0); 
    for tl_rxn_idx = 1:length(tl_elong_rxn_name_list)
        tl_rxn_obj = sbioselect(Mobj,'Type','Reaction','Reaction',tl_elong_rxn_name_list{tl_rxn_idx}); 
        for tl_rxn_obj_idx = 1:length(tl_rxn_obj)
            tl_rxn_obj_ind = tl_rxn_obj(tl_rxn_obj_idx); 
            
            addproduct(tl_rxn_obj_ind,'toxin',1); 
        end
    end
        % Add a reaction for toxin buffering
    toxin_removal = addreaction(Mobj,'toxin -> null'); 
    set(toxin_removal, 'ReactionRate', 'tx_capacity_param');
        % Add a reaction for Ribo degradation 
    Ribo_deg = addreaction(Mobj,'Ribo -> null');
    set(Ribo_deg,'ReactionRate','k_toxin * (1 + tanh(toxin - toxin_threshold)) * Ribo')



    % Export Mobj.Parameters to an excel file if there's not one.
    % Otherwise load parameter info from the excel file 
 
    if ~exist(output_path,'file')
        save_parameter_info(Mobj,output_path); 
    end
    parameter_info_table = readtable(output_path); 
    % check existence of promotor, check existence of empty gene expression 
        % Initialize everything to be not present 
    T7_strong_flag = false;
    T7_weak_flag = false; 
    sigma70_strong_flag = false;
    sigma70_weak_flag = false; 
    empty_flag = false; 
    
    for description_idx = 1:length(group_description)
        group_description_single = group_description{description_idx}; 
        if ~contains(group_description,'no_empty')
            empty_flag = true; 
        end
        if contains(group_description_single,'sigma70_weak')
            sigma70_weak_flag = true;
        elseif contains(group_description_single,'sigma70')
            sigma70_strong_flag = true;
        end
        if contains(group_description_single,'T7_weak')
            T7_weak_flag = true;
        elseif contains(group_description_single,'T7')
            T7_strong_flag = true;
        end
    end

    % Add parameters to be estimated in estimatedInfo
        % Compile a cell of parameter names
        paramsToEstimate = {};
        initialValues = []; 
        bounds = []; % n * 2 matrix 
        model_rules = get(Mobj,'Rules'); 
        % loop through Mobj.Parameters and append names 
        Mobj_parameters = get(Mobj,'Parameters'); 
        for param_idx = 1:length(Mobj_parameters)
            Mobj_param = Mobj_parameters(param_idx);
            add_flag = true; % Whether the parameter should be estimated, default to be true 
    
            % Don't estimate the '_F' parameters and just fix them as 1's 
            if contains(Mobj_param.Name,'_F')
                add_flag = false;
            end

            % Don't estimate AGTPreg_varying 
            if contains(Mobj_param.Name,'AGTPreg_varying')
                add_flag = false;
                ruleObj = addrule(Mobj,'AGTPreg_varying=AGTPreg_ON'); % add initial assignment of AGTP regeneration rate
            end
    
            % Revise the rule assignment on AGTP
            % if strcmp(Mobj_param.Name,'AGTPdeg_time')
            %     add_flag = false;
            % end

            % If parameter is assigned in Model.Rules, no need to estimate 
            for rule_idx = 1:length(model_rules)
                rule = model_rules(rule_idx).Rule; 
                if contains(rule(1:find(rule=='=')),Mobj_param.Name) && ~contains(Mobj_param.Name,'_0')
                    add_flag = false; 
                end
            end
            
            % If parameter is not relevant to data, don't add it
            if contains(Mobj_param.Name,'t7') && ~T7_weak_flag && ~T7_strong_flag
                add_flag = false;
            end
            if contains(Mobj_param.Name,'T7') && ~T7_weak_flag && ~T7_strong_flag
                add_flag = false;
            end
            if contains(Mobj_param.Name,'PT773') && ~T7_weak_flag
                add_flag = false; 
            end
            if contains(Mobj_param.Name,'PT7') && ~T7_strong_flag
                add_flag = false; 
            end
            if contains(Mobj_param.Name,'J23119') && ~sigma70_strong_flag
                add_flag = false;
            end
            if contains(Mobj_param.Name,'J23105') && ~sigma70_weak_flag
                add_flag = false;
            end

            % If there's only PE data, don't estimate broc-related
            % parameters
            if isempty(tx_selected_table) && contains(Mobj_param.Name,'broc')
                add_flag = false; 
            end
    
            table_idx = find(strcmp(parameter_info_table.Name,Mobj_param.Name));
            LB = table2array(parameter_info_table(table_idx,'LB'));
            UB = table2array(parameter_info_table(table_idx,'UB'));
            %If parameter is fixed (i.e. LB=UB), no need to estimate 
            if isequal(LB,UB)
                add_flag = false;
            end
            if add_flag
                paramsToEstimate{end + 1} = strcat('log(',Mobj_param.Name,')');
                initVal = table2array(parameter_info_table(table_idx,'InitVal')); 
                initialValues = [initialValues initVal];
                bounds = [bounds;[LB,UB]]; 
                
            end
        end

    % Update model in problemObject
    problemObject.Model = Mobj; 
    problemObject.Estimated = estimatedInfo(paramsToEstimate,'InitialValue',initialValues,'Bounds',bounds);
    % problemObject.FunctionName = 'fmincon'; 
    problemObject.Pooled = true; 
    problemObject.ProgressPlot = 0; 


%% Count number of parameters in the model 
fprintf('Number of parameters to be estimated:%d',length(paramsToEstimate))

end