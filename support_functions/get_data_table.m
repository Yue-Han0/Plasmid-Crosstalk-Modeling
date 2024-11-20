function grouped_data = get_data_table(group_description)
    %% Load data 
        % load both transcription and protein expression data 
    tx_data_file = load('data_structures/simbio_data_table_updated_FP.mat'); 
    PE_data_file = load('data_structures/simbio_data_table_PE_updated_FP.mat');
    tx_data_table = tx_data_file.all_data_table; 
    PE_data_table = PE_data_file.all_data_table; 
    tx_data_description = tx_data_file.variable_name_stem_list; 
    PE_data_description = PE_data_file.variable_name_stem_list; 

        % Select relevant data 
    [tx_group_number,tx_selected_table] = get_relevant_data(tx_data_table,tx_data_description,group_description); 
    [PE_group_number,PE_selected_table] = get_relevant_data(PE_data_table,PE_data_description,group_description); 

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
    grouped_data = groupedData(selected_table); % define a groupedData object 
    grouped_data.Properties.IndependentVariableName = 'Time'; 
    grouped_data.Properties.GroupVariableName = 'Group'; 

end
