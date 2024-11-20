clear
clc
% 03/10/2023 edit: add a column of standard deviation next to fluorescence
% for confidence interval calculation 

% 03/29/2023 edit: add a column with gain 

% 06/12/2023 edit: generate data table for protein data 

% 10/10/2023 edit: Change empty plasmid concentration 
%% Definitions

% promotor_name_list = {'T7_strong','T7_weak','sigma70'};
% for the updated data 
% promotor_name_list = {'T7_strong','T7_weak','sigma70_strong'};
promotor_name_list = {'T7_strong','T7_weak','sigma70_strong','sigma70_weak'}; % protein data 
mix_name_list = {'no_empty','empty','empty_T7','empty_sigma70'};
appendix1 = '_updated_FP'; 
appendix2 = '_updated_AP'; 
% appendix2 = ''; 
sigma_appendix_list = {'','_lowConc','_highConc'};
reporter_conc_list = [0.5,1,2.5,5,10,15,30]; 
reporter_lowConc_list = reporter_conc_list(1:4);
reporter_highConc_list = reporter_conc_list(5:7);
reporter_conc_index_sigma_appendix = {reporter_conc_list,reporter_lowConc_list,reporter_highConc_list}; 
gain_index_sigma_appendix = [70,100,100];

%% Load data 
% % T7 strong
%     % No empty 
% T7_strong_no_empty_data = readtable('processed_data/T7_strong_no_empty.xlsx',Sheet='fluorescence_data');
% T7_strong_no_empty_timeVec = readtable('processed_data/T7_strong_no_empty.xlsx',Sheet='timeVec');
% T7_strong_no_empty_stdev = readtable('processed_data/T7_strong_no_empty.xlsx',Sheet='stdev');
%     % empty 
% T7_strong_empty_data = readtable('processed_data/T7_strong_empty.xlsx',Sheet='fluorescence_data');
% T7_strong_empty_timeVec = readtable('processed_data/T7_strong_empty.xlsx',Sheet='timeVec');
% T7_strong_empty_stdev = readtable('processed_data/T7_strong_empty.xlsx',Sheet='stdev');
%     % empty T7
% T7_strong_empty_T7_data = readtable('processed_data/T7_strong_empty_T7.xlsx',Sheet='fluorescence_data');
% T7_strong_empty_T7_timeVec = readtable('processed_data/T7_strong_empty_T7.xlsx',Sheet='timeVec');
% T7_strong_empty_T7_stdev = readtable('processed_data/T7_strong_empty_T7.xlsx',Sheet='stdev');
%     % empty sigma70
% T7_strong_empty_sigma70_data = readtable('processed_data/T7_strong_empty_sigma70.xlsx',Sheet='fluorescence_data');
% T7_strong_empty_sigma70_timeVec = readtable('processed_data/T7_strong_empty_sigma70.xlsx',Sheet='timeVec');
% T7_strong_empty_sigma70_stdev = readtable('processed_data/T7_strong_empty_sigma70.xlsx',Sheet='stdev');
% %T7 weak
%     % No empty 
% T7_weak_no_empty_data = readtable('processed_data/T7_weak_no_empty.xlsx',Sheet='fluorescence_data');
% T7_weak_no_empty_timeVec = readtable('processed_data/T7_weak_no_empty.xlsx',Sheet='timeVec');
% T7_weak_no_empty_stdev = readtable('processed_data/T7_weak_no_empty.xlsx',Sheet='stdev');
%     % empty
% T7_weak_empty_data = readtable('processed_data/T7_weak_empty.xlsx',Sheet='fluorescence_data');
% T7_weak_empty_timeVec = readtable('processed_data/T7_weak_empty.xlsx',Sheet='timeVec');
% T7_weak_empty_stdev = readtable('processed_data/T7_weak_empty.xlsx',Sheet='stdev');
%     % empty T7
% T7_weak_empty_T7_data = readtable('processed_data/T7_weak_empty_T7.xlsx',Sheet='fluorescence_data');
% T7_weak_empty_T7_timeVec = readtable('processed_data/T7_weak_empty_T7.xlsx',Sheet='timeVec');
% T7_weak_empty_T7_stdev = readtable('processed_data/T7_weak_empty_T7.xlsx',Sheet='stdev');
%     % empty sigma70 
% T7_weak_empty_sigma70_data = readtable('processed_data/T7_weak_empty_sigma70.xlsx',Sheet='fluorescence_data');
% T7_weak_empty_sigma70_timeVec = readtable('processed_data/T7_weak_empty_sigma70.xlsx',Sheet='timeVec');
% T7_weak_empty_sigma70_stdev = readtable('processed_data/T7_weak_empty_sigma70.xlsx',Sheet='stdev');
% % sigma70 strong
%     % Low conc
%     % No empty 
% sigma70_strong_no_empty_data_lowConc = readtable('processed_data/sigma70_no_empty_lowConc.xlsx',Sheet='fluorescence_data');
% sigma70_strong_no_empty_timeVec_lowConc = readtable('processed_data/sigma70_no_empty_lowConc.xlsx',Sheet='timeVec');
% sigma70_strong_no_empty_stdev_lowConc = readtable('processed_data/sigma70_no_empty_lowConc.xlsx',Sheet='stdev');
%     % empty 
% sigma70_strong_empty_data_lowConc = readtable('processed_data/sigma70_empty_lowConc.xlsx',Sheet='fluorescence_data');
% sigma70_strong_empty_timeVec_lowConc = readtable('processed_data/sigma70_empty_lowConc.xlsx',Sheet='timeVec');
% sigma70_strong_empty_stdev_lowConc = readtable('processed_data/sigma70_empty_lowConc.xlsx',Sheet='stdev');
%     % empty T7
% sigma70_strong_empty_T7_data_lowConc = readtable('processed_data/sigma70_empty_T7_lowConc.xlsx',Sheet='fluorescence_data');
% sigma70_strong_empty_T7_timeVec_lowConc = readtable('processed_data/sigma70_empty_T7_lowConc.xlsx',Sheet='timeVec');
% sigma70_strong_empty_T7_stdev_lowConc = readtable('processed_data/sigma70_empty_T7_lowConc.xlsx',Sheet='stdev');
%     % empty sigma70
% sigma70_strong_empty_sigma70_data_lowConc = readtable('processed_data/sigma70_empty_sigma70_lowConc.xlsx',Sheet='fluorescence_data');
% sigma70_strong_empty_sigma70_timeVec_lowConc = readtable('processed_data/sigma70_empty_sigma70_lowConc.xlsx',Sheet='timeVec');
% sigma70_strong_empty_sigma70_stdev_lowConc = readtable('processed_data/sigma70_empty_sigma70_lowConc.xlsx',Sheet='stdev');
%     % High conc
%     % No empty 
% sigma70_strong_no_empty_data_highConc = readtable('processed_data/sigma70_no_empty_highConc.xlsx',Sheet='fluorescence_data');
% sigma70_strong_no_empty_timeVec_highConc = readtable('processed_data/sigma70_no_empty_highConc.xlsx',Sheet='timeVec');
% sigma70_strong_no_empty_stdev_highConc = readtable('processed_data/sigma70_no_empty_highConc.xlsx',Sheet='stdev');
%     % empty 
% sigma70_strong_empty_data_highConc = readtable('processed_data/sigma70_empty_highConc.xlsx',Sheet='fluorescence_data');
% sigma70_strong_empty_timeVec_highConc = readtable('processed_data/sigma70_empty_highConc.xlsx',Sheet='timeVec');
% sigma70_strong_empty_stdev_highConc = readtable('processed_data/sigma70_empty_highConc.xlsx',Sheet='stdev');
%     % empty T7
% sigma70_strong_empty_T7_data_highConc = readtable('processed_data/sigma70_empty_T7_highConc.xlsx',Sheet='fluorescence_data');
% sigma70_strong_empty_T7_timeVec_highConc = readtable('processed_data/sigma70_empty_T7_highConc.xlsx',Sheet='timeVec');
% sigma70_strong_empty_T7_stdev_highConc = readtable('processed_data/sigma70_empty_T7_highConc.xlsx',Sheet='stdev');
%     % empty sigma70
% sigma70_strong_empty_sigma70_data_highConc = readtable('processed_data/sigma70_empty_sigma70_highConc.xlsx',Sheet='fluorescence_data');
% sigma70_strong_empty_sigma70_timeVec_highConc = readtable('processed_data/sigma70_empty_sigma70_highConc.xlsx',Sheet='timeVec');
% sigma70_strong_empty_sigma70_stdev_highConc = readtable('processed_data/sigma70_empty_sigma70_highConc.xlsx',Sheet='stdev');
% 


% A version without these hard-code loading (using eval and statements in strings) 
dataset_idx = 1; 
variable_name_stem_list = {};

for promotor_idx = 1:length(promotor_name_list)
    promotor_name = promotor_name_list{promotor_idx}; 
    for mix_idx = 1:length(mix_name_list)
        mix_name = mix_name_list{mix_idx}; 
        % define concentration of DNAs' in system 
        if contains(mix_name,'empty_T7') || contains(mix_name,'empty_sigma70')
            empty_conc = 10; 
        else
            empty_conc = 0; 
        end
        if ~contains(mix_name,'no_empty')
            kan_conc = 10;
        else
            kan_conc = 0; 
        end

%         % Old data
%         if contains(promotor_name,'sigma70')
%             sigma_appendix_loop = 2:3;
%         else
%             sigma_appendix_loop = 1;
%         end
        % updated data 
        sigma_appendix_loop = 1; 

        for sigma_appendix_idx = sigma_appendix_loop
            sigma_appendix = sigma_appendix_list{sigma_appendix_idx};
            plasmid_doses = reporter_conc_index_sigma_appendix{sigma_appendix_idx}; 
            % load data into workspace 

                % Old data 
%             [statement_data,statement_timeVec,statement_stdev,variable_name_stem]= get_statements(promotor_name,mix_name,sigma_appendix); 
                % Updated data 
            [statement_data,statement_timeVec,statement_stdev,variable_name_stem]= get_statements(promotor_name,mix_name,sigma_appendix,appendix2); 

            eval(statement_data)
            eval(statement_timeVec)
            eval(statement_stdev)
            % organize variable names 
            [variable_name_stem_list{end + 1:end + length(plasmid_doses),1}] = deal(variable_name_stem);
            % construct data table with currently loaded data 
            construct_data_table_statement = sprintf('[data_table,dataset_idx] = construct_data_table(%s_data,%s_stdev,%s_timeVec,plasmid_doses,dataset_idx,kan_conc,empty_conc,promotor_name,mix_name);',...
                variable_name_stem,variable_name_stem,variable_name_stem); 
            eval(construct_data_table_statement)
            if ~exist('all_data_table','var')
                all_data_table = data_table;
            else
                all_data_table = [all_data_table;data_table];
            end
        end
    end
end

% Let's just calculate CI instead of trying to use normal distribution
% objects
% CI = mean +- z * (sd/sqrt(N))
fluorescence_data = all_data_table.fluorescence;
stdev_data = all_data_table.stdev; 
N = 3; 
z_score = 1.96; % for 95% CI
CI_lb = fluorescence_data - z_score * (stdev_data/sqrt(N));
CI_ub = fluorescence_data + z_score * (stdev_data/sqrt(N));
all_data_table.CI_lb = CI_lb; 
all_data_table.CI_ub = CI_ub; 

% Add gain in the data table 
%     % AP lysate 
%     % Gain = 70 if sigma70 promotor; gain = 75 if T7 promotor 
%     % Define the group number that separates T7 and sigma70
% sigma70_starting_group = 57; 
% starting_group_idx_list_datatable = find(all_data_table.Group==sigma70_starting_group); 
% starting_group_idx_datatable = starting_group_idx_list_datatable(1); 
% gain = nan(height(all_data_table),1);
% gain(1:starting_group_idx_datatable-1) = repmat(75,starting_group_idx_datatable-1,1);
% gain(starting_group_idx_datatable:end) = repmat(70,length(gain)-starting_group_idx_datatable+1,1); 
    % FP lysate 
    % Gain = 70 if sigma70 promotor & T773; gain = 75 if T7 strong promotor
    % Define the group number that separates T7 strong and T7 weak 
% T7weak_starting_group = 29; 
% starting_group_idx_list_datatable = find(all_data_table.Group==T7weak_starting_group); 
% starting_group_idx_datatable = starting_group_idx_list_datatable(1); 
% gain = nan(height(all_data_table),1);
% gain(1:starting_group_idx_datatable-1) = repmat(75,starting_group_idx_datatable-1,1);
% gain(starting_group_idx_datatable:end) = repmat(70,length(gain)-starting_group_idx_datatable+1,1); 

% Add PE fluorescence data gain in data table 
sigma70_weak_starting_group = 85; 
starting_group_idx_list_datatable = find(all_data_table.Group==sigma70_weak_starting_group); 
starting_group_idx_datatable = starting_group_idx_list_datatable(1); 
gain = nan(height(all_data_table),1);
gain(1:starting_group_idx_datatable-1) = repmat(50,starting_group_idx_datatable-1,1);
gain(starting_group_idx_datatable:end) = repmat(75,length(gain)-starting_group_idx_datatable+1,1); 

% Add var to table 
all_data_table = addvars(all_data_table,gain);
save('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/benchmark_model/data_structures/simbio_data_table_PE_updated_FP','all_data_table','variable_name_stem_list'); 

function [statement_data,statement_timeVec,statement_stdev,variable_name_stem] = get_statements(promotor_name,mix_name,sigma_appendix,appendix)
    if ~exist('appendix','var')
        appendix = ''; 
    end
    statement_data = sprintf('PE_%s_%s%s_data = readtable(''processed_data/PE_%s_%s%s%s.xlsx'',''ReadVariableNames'',true,Sheet=''fluorescence_data'');',promotor_name,mix_name,sigma_appendix,promotor_name,mix_name,sigma_appendix,appendix);
    statement_timeVec = sprintf('PE_%s_%s%s_timeVec = readtable(''processed_data/PE_%s_%s%s%s.xlsx'',''ReadVariableNames'',true,Sheet=''timeVec'');',promotor_name,mix_name,sigma_appendix,promotor_name,mix_name,sigma_appendix,appendix);
    statement_stdev = sprintf('PE_%s_%s%s_stdev = readtable(''processed_data/PE_%s_%s%s%s.xlsx'',''ReadVariableNames'',true,Sheet=''stdev'');',promotor_name,mix_name,sigma_appendix,promotor_name,mix_name,sigma_appendix,appendix);

    variable_name_stem = sprintf('PE_%s_%s%s',promotor_name,mix_name,sigma_appendix);
end

function [data_table,dataset_idx] = construct_data_table(data,stdev,timeVec,plasmid_doses,dataset_idx,kan_conc,empty_conc,promotor_name,mix_name)
    %% Build a table with Time, 3 kind of plasmid concentrations, and fluorescence, group name, group idx

    % sanity check
    if ~isequal(size(data,2) - 1, length(plasmid_doses))
        error('Inconsistent dose & response size')
    end

    % Add timeVec to data table and remove index column 
    data_table = addvars(data,timeVec.Time,'NewVariableNames',{'Time'},'Before',1);
    data_table = removevars(data_table,'Var1');

    % Let's organize the data in a different way 
        % Let's have 9 columns: timeVec, fluorescence, groupName (aka dataset_idx),
        % T7_strong_reporter_plasmid, T7_weak_reporter_plasmid,
        % sigma70_reporter_plasmid, kan_plasmid, empty_T7_plasmid,
        % empty_sigma70 plasmid
    augmented_timeVec = repmat(timeVec.Time,[length(plasmid_doses),1]);

    % Augmented others should have the same dimension as augmented_timeVec
    augmented_fluorescence = nan(length(augmented_timeVec),1);
    augmented_stdev = nan(length(augmented_timeVec),1); 
    groupName = nan(length(augmented_timeVec),1);
    reporter_plasmid = zeros(length(augmented_timeVec),1); % a 'glossary' term for all reporter plasmids 
    T7_strong_reporter_plasmid = zeros(length(augmented_timeVec),1); 
    T7_weak_reporter_plasmid = zeros(length(augmented_timeVec),1); 
    sigma70_strong_reporter_plasmid = zeros(length(augmented_timeVec),1); 
    sigma70_weak_reporter_plasmid = zeros(length(augmented_timeVec),1); 
    kan_plasmid = zeros(length(augmented_timeVec),1); 
    empty_plasmid = zeros(length(augmented_timeVec),1); % glossary term for empty plasmids 
    empty_T7_plasmid = zeros(length(augmented_timeVec),1); 
    empty_sigma70_plasmid = zeros(length(augmented_timeVec),1); 
    
    variable_names = data.Properties.VariableNames; 
    for conc_idx = 1:length(plasmid_doses)
        reporter_dna_conc = plasmid_doses(conc_idx); 
        kan_dna_conc = kan_conc + reporter_dna_conc; 
        fluorescence_data = data{:,variable_names{conc_idx+1}}; 
        stdev_data = stdev{:,variable_names{conc_idx+1}};

        augmented_fluorescence((conc_idx - 1) * length(fluorescence_data)+1:conc_idx * length(fluorescence_data),1) = ...
            fluorescence_data;
        augmented_stdev((conc_idx - 1) * length(stdev_data)+1:conc_idx * length(stdev_data),1) = ...
            stdev_data;
        groupName((conc_idx - 1) * length(fluorescence_data)+1:conc_idx * length(fluorescence_data)) = ...
            repmat(dataset_idx+conc_idx-1,[length(fluorescence_data),1]); 
        reporter_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = reporter_dna_conc;
        if strcmp(promotor_name,'T7_strong')
            T7_strong_reporter_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = reporter_dna_conc;
        elseif strcmp(promotor_name,'T7_weak')
            T7_weak_reporter_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = reporter_dna_conc;
        elseif strcmp(promotor_name,'sigma70_strong')
            sigma70_strong_reporter_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = reporter_dna_conc;
        elseif strcmp(promotor_name,'sigma70_weak')
            sigma70_weak_reporter_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = reporter_dna_conc;
        end
        kan_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = kan_dna_conc;
        empty_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = empty_conc;
        if strcmp(mix_name,'empty_T7')
            empty_T7_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = empty_conc;
        elseif strcmp(mix_name,'empty_sigma70')
            empty_sigma70_plasmid((conc_idx - 1) * length(fluorescence_data) + 1) = empty_conc;
        end

    end
    
    data_table = table(groupName,T7_strong_reporter_plasmid,T7_weak_reporter_plasmid,sigma70_strong_reporter_plasmid,...
        sigma70_weak_reporter_plasmid,kan_plasmid,empty_T7_plasmid,empty_sigma70_plasmid,reporter_plasmid,...
        empty_plasmid,augmented_timeVec,augmented_fluorescence,augmented_stdev,...
        'VariableNames',{'Group','T7_strong_GFP_plasmid','T7_weak_GFP_plasmid',...
        'sigma70_strong_GFP_plasmid','sigma70_weak_GFP_plasmid','kan_plasmid','empty_T7_plasmid','empty_sigma70_plasmid',...
        'GFP_plasmid','empty_plasmid','Time','fluorescence','stdev'});
    dataset_idx = dataset_idx + length(plasmid_doses); 

end







