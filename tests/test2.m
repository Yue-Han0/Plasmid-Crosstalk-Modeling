% clear
% clc


currentpath = pwd;
addpath(genpath(sprintf('%s/plasmid_crosstalk_config_files',currentpath))); 
addpath(genpath(sprintf('%s/PESTO',currentpath))); 
addpath(genpath(sprintf('%s/txtlsim_buildacell',currentpath))); 
addpath(genpath(sprintf('%s',currentpath)));
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath));
% 
promotors = {'T7 strong','T7 weak','sigma70 strong','sigma70 weak'};
unique_empty_description = {'no empty','empty','empty T7','empty sigma70'};
color_list = {'g',[0.5,0.5,0.5],'r','b'};
conc_vec = [0.5,1,2.5,5,10,15,30]; 
num_conc = length(conc_vec); 
% num_promotor = length(promotors); 
% empty_conc = 15; 

%% Check missing large-scale fitting result file
% timestamp_confusion = {'1531','1532','1533','1534'}; 
% result_fileName_stem = 'param_est_run_save/20231102_param_est_run';
% 
% for run_idx = 1:100
%     result_fileName_suffix = sprintf('_%d.mat',run_idx); 
%     file_exist_flag = 0; 
%     for timeStamp_idx = 1:length(timestamp_confusion)
%         timestamp = timestamp_confusion{timeStamp_idx}; 
%         test_result_fileName = strcat(result_fileName_stem,timestamp,result_fileName_suffix); 
%         if exist(test_result_fileName,'file')
%             file_exist_flag = 1;
%             break
%         end
%     end
% 
%     if ~file_exist_flag
%         disp(run_idx)
%         run_number = run_idx + 100; 
%         % If doesn't exist, generate a job submission file with new
%         % initial seeds and add to a master job submission txt file
%         command = fopen(sprintf('jobSubmission/more_initial_points_%d.txt',run_number),'w');
%         fprintf(command,'#!/bin/bash \n #SBATCH -JmoreInitialPoint \n #SBATCH --account=gts-mstyczynski6 \n #SBATCH -N1 -n1 \n #SBATCH -t5760 \n #SBATCH -qinferno #SBATCH -oReport-%%j.out \n cd $SLURM_SUBMIT_DIR \n module load matlab');
%         fprintf(command,'\n srun matlab -nosplash -nodisplay -singleCompThread -r "wrapper_fit_to_data_PACE(%d)"',run_number);
%         fclose('all');
% 
%             % fprintf(command,'\n sbatch -A gts-mstyczynski6 -q inferno -N1 --ntasks-per-node=1 --mem-per-cpu=5G -t 96:00:00 jobSubmission/more_initial_points_%d.txt',run_number);
%     end
% end

%% 11/1 Test toxin mechanism 
    % Result file 
% ctrl_TX_gen_on_RNAP_resultFileName ='param_est_run_save/20231101_param_est_run0929_RNAP_toxin_mech.mat';
% TX_gen_on_Ribo_resultFileName = 'param_est_run_save/20230125_param_est_run1111_TXgen_Ribo.mat';
% TL_gen_on_RNAP_resultFileName = 'param_est_run_save/20230125_param_est_run2230_TLgen_RNAP.mat';
% TL_gen_on_Ribo_resultFileName = 'param_est_run_save/20230125_param_est_run1326_TLgen_Ribo.mat';
% 
% load(TX_gen_on_Ribo_resultFileName); 
% % load param_est_run_save/20230813_param_est_run0807_PE_T7_strong_no_empty.mat
% % Get the SSEs 
% all_SSE = nan(length(all_fitResults),1); 
% for i = 1:length(all_fitResults)
%    fitResult = all_fitResults{i,1}; 
%    all_SSE(i) = fitResult.SSE;
% end
% 
% [~,opt_idx] = min(all_SSE); 
% opt_fitResult = all_fitResults{opt_idx,1}; 
% % 
% fitted_data = fitted(opt_fitResult); 
% 
% [Time,Data] = process_crosstalk_data_from_source(fitted_data,'simulated_data');
% % [Time,Data] = process_crosstalk_data_from_source(problemObject.Data,'grouped_data');
% figure; 
% for conc_idx = 1:length(Data)
%     plot(Time{conc_idx}./3600,Data{conc_idx},'LineWidth',1.5)
%     hold on 
% end
% sgtitle('TL-gen toxin on Ribosome')
% xlabel('Time (hr)','FontSize',14)
% ylabel('Concentration (nM)','FontSize',14)
% legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')
% set(gca,'FontSize',14)

%% 10/30 Test toxin mechanism after placing toxin generation on consumption rather than elongation 
% load param_est_run_save/20231031_param_est_run1051_RNAP_toxin_mech.mat
    % Plot out all 48 fitting results 
% for iter = 1:length(all_fitResults)
%     fitResult = all_fitResults{iter,1};
%     simData = fitted(fitResult); 
%     [Time,Data] = process_crosstalk_data_from_source(simData,'simulated_data');
%     subplot(7,7,iter)
%     for conc_idx = 1:length(Data)
%         plot(Time{conc_idx},Data{conc_idx},'LineWidth',1.5)
%         hold on 
%     end
%     title(sprintf('iter #%d',iter))
% end
% sgtitle('RNAP toxin mechanism fitting')
% % saveas(gcf,'plots/20231030_RNAP_toxin_mechanism_all_fitting','png'); 

% Get the SSEs 
% all_SSE = nan(length(all_fitResults),1); 
% for i = 1:length(all_fitResults)
%    fitResult = all_fitResults{i,1}; 
%    all_SSE(i) = fitResult.SSE;
% end
% 
% [~,opt_idx] = min(all_SSE); 
% opt_fitResult = all_fitResults{opt_idx,1}; 
% % 
% fitted_data = fitted(opt_fitResult); 

%% Check t_broc_max in different conditions 
%     % Load tx data 
% tx_data_file = load('data_structures/simbio_data_table_updated_FP.mat'); 
% tx_data_table = tx_data_file.all_data_table; 
% tx_data_description = tx_data_file.variable_name_stem_list; 
% group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% [tx_group_number,tx_selected_table] = get_relevant_data(tx_data_table,tx_data_description,group_description);
% Time = [tx_selected_table.Time];
% fluorescence = [tx_selected_table.fluorescence];
% time_start_idx_list = find(Time < 90);
% max_time_list = nan(84,1); 
% max_diff_time_list = nan(84,1); 
% for j = 1:length(time_start_idx_list)
%     time_start_idx = time_start_idx_list(j);
%     if isequal(j,length(time_start_idx_list))
%         time_end_idx = length(Time);
%     else
%         time_end_idx = time_start_idx_list(j + 1) - 1; 
%     end
%     single_timeVec = Time(time_start_idx + 1:time_end_idx); 
%     single_fluorescence_time_course = fluorescence(time_start_idx + 1 :time_end_idx); % plus 1 to avoid outlier at early timepoints 
%     diff_single_fluorescence = diff(single_fluorescence_time_course); 
%     [~,max_increase_idx] = max(diff_single_fluorescence); 
%     [~,max_idx] = max(single_fluorescence_time_course); 
%     max_time = single_timeVec(max_idx);
%     max_time_list(j) = max_time;
%     max_diff_time = single_timeVec(max_increase_idx);
%     max_diff_time_list(j) = max_diff_time; 
% end
% 
%     % Plot max time as a function of increasing reporter plasmid concentration 
% for prom_idx = 1:length(promotors) - 1
%     subplot(2,2,prom_idx)
%     title(promotors{prom_idx})
%     for empty_idx = 1:length(unique_empty_description)
%         start_idx = (prom_idx - 1) * length(unique_empty_description) * length(conc_vec) + ...
%             (empty_idx - 1) * length(conc_vec) + 1; 
%         end_idx = (prom_idx - 1) * length(unique_empty_description) * length(conc_vec) + ...
%             empty_idx * length(conc_vec); 
%         color = color_list{empty_idx}; 
%         plot(conc_vec,max_diff_time_list(start_idx:end_idx),'LineWidth',1.5,'Color',color); 
%         hold on 
%         xlabel('Time(s)')
%         ylabel('Peak time')
%     end
% end
% legend('no empty','empty','empty T7','empty sigma70')
% sgtitle('Max Broc increase time')

%% Build a Simbiology model for Broc mRNA production and degradation & perform parameter estimation 
%     % Load data for mRNA production 
% tx_data_file = load('data_structures/simbio_data_table_updated_FP.mat'); 
% tx_data_table = tx_data_file.all_data_table; 
% tx_data_description = tx_data_file.variable_name_stem_list; 
% group_description = {'sigma70_strong_no_empty'}; 
%     % Select relevant data 
% [tx_group_number,tx_selected_table] = get_relevant_data(tx_data_table,tx_data_description,group_description); 
% 
%     % Convert fluorescence to concentration and add concentration to
%     % data table 
% if ~isempty(tx_selected_table)
%     split_point_idx_list = find(tx_selected_table.gain==tx_selected_table.gain(1)); 
%     split_point_idx = split_point_idx_list(end); 
%     if isequal(length(split_point_idx_list),height(tx_selected_table))% if there's only one gain 
%         mRNA_concentration = [fluo_mRNA_conversion(tx_selected_table.fluorescence,[],tx_selected_table.gain(1))];
%     else
%             % Applies to when there are 2 gains 
%         mRNA_concentration = [fluo_mRNA_conversion(tx_selected_table.fluorescence(1:split_point_idx),[],tx_selected_table.gain(1));...
%             fluo_mRNA_conversion(tx_selected_table.fluorescence(split_point_idx+1:end),[],tx_selected_table.gain(split_point_idx + 1))];
%     end
%     tx_selected_table = addvars(tx_selected_table,mRNA_concentration,'NewVariableNames','mRNA_concentration'); 
%         % Replace negative values with nan 
%     tx_selected_table = remove_neg_val(tx_selected_table);
%         % Find new split after negative value removal
%     new_split_point_idx_list = find(tx_selected_table.gain==tx_selected_table.gain(1)); 
%     new_split_point_idx = new_split_point_idx_list(end); 
%     if isequal(length(new_split_point_idx_list),height(tx_selected_table))
%         mRNA_CI_lb = [fluo_mRNA_conversion(tx_selected_table.CI_lb,[],tx_selected_table.gain(1))];
%         mRNA_CI_ub = [fluo_mRNA_conversion(tx_selected_table.CI_ub,[],tx_selected_table.gain(1))];
%     else
%         mRNA_CI_lb = [fluo_mRNA_conversion(tx_selected_table.CI_lb(1:new_split_point_idx),[],tx_selected_table.gain(1));...
%             fluo_mRNA_conversion(tx_selected_table.CI_lb(new_split_point_idx+1:end),[],tx_selected_table.gain(new_split_point_idx + 1))];
%         mRNA_CI_ub = [fluo_mRNA_conversion(tx_selected_table.CI_ub(1:new_split_point_idx),[],tx_selected_table.gain(1));...
%             fluo_mRNA_conversion(tx_selected_table.CI_ub(new_split_point_idx+1:end),[],tx_selected_table.gain(new_split_point_idx + 1))];
%     end
%     tx_selected_table = addvars(tx_selected_table,mRNA_CI_lb,'NewVariableNames','mRNA_concentration_CI_lb'); 
%     tx_selected_table = addvars(tx_selected_table,mRNA_CI_ub,'NewVariableNames','mRNA_concentration_CI_ub'); 
% end
% 
% 
% % Add to problemObject
% RNA_prod_deg_probObject = fitproblem(); % Initialize a fitProblem object
%         % Data
% RNA_prod_deg_probObject.Data = groupedData(tx_selected_table); % define a groupedData object 
% RNA_prod_deg_probObject.Data.Properties.IndependentVariableName = 'Time'; 
% RNA_prod_deg_probObject.Data.Properties.GroupVariableName = 'Group'; 
%         % Model 
% Model_struct = sbioloadproject('simbiology_models/mRNA_prod_deg.sbproj'); 
% RNA_prod_deg_Model = Model_struct.m1; 
%             % Process the model to make sure all reaction kinetic
%             % parameters are added to model and deleted in rxn kinetic law 
% Mobj_rxns = get(RNA_prod_deg_Model,'Reactions');
% for rxn_idx = 1:length(Mobj_rxns)
%     rxn = Mobj_rxns(rxn_idx);
%     rxn_params = rxn.KineticLaw.Parameters;
%     for rxn_param_idx = 1:length(rxn_params)
%         rxn_param = rxn_params(rxn_param_idx);
%         rxn_param_name = rxn_param.Name; 
%         % Is this parameter already in model parameters?
%         global_param = sbioselect(RNA_prod_deg_Model.Parameters,'Name',rxn_param_name);
%         if ~isempty(global_param)
%             % if already in model parameters, select and delete
%             fprintf('\n %s already exists in model object,deleting',rxn_param_name)
%             delete(rxn_param); 
%         else
%             % If doesn't exist in model parameters, move to model
%             % parameters 
%                 % Make sure parameter name matches 
%             if ~strcmp(rxn_param_name,rxn_param.Name)
%                 %If doesn't match, set parameter name
%                 set(rxn_param,'Name',rxn_param_name); 
%             end
%             test_obj = move(rxn_param,RNA_prod_deg_Model);
%         end
%     end
% end
%             % Set initial condition for AGTP & CUTP
% AGTP_species = sbioselect(RNA_prod_deg_Model.Species,'Name','AGTP'); 
% set(AGTP_species,'Value',2.05e+06); 
% CUTP_species = sbioselect(RNA_prod_deg_Model.Species,'Name','CUTP'); 
% set(CUTP_species,'Value',1.7e+06); 
% 
% RNA_prod_deg_probObject.Model = RNA_prod_deg_Model; 
% RNA_prod_deg_probObject.ResponseMap = '[RNA utrbroc--no_protein] = mRNA_concentration';
% 
% 
%         % Doses
% RNA_prod_deg_probObject.Doses = createDoses(RNA_prod_deg_probObject.Data,{'sigma70_strong_reporter_plasmid'}); 
% for dose_idx = 1:size(RNA_prod_deg_probObject.Doses,1)
%     RNA_prod_deg_probObject.Doses(dose_idx,1).TargetName = 'DNA pJ23119--utrbroc--no_protein';
% end
%         % Estimates 
% param_info_table = readtable('parameters_mRNA_prod_deg.xlsx');
% Mobj_parameters = get(RNA_prod_deg_Model,'Parameters'); 
% model_rules = get(RNA_prod_deg_Model,'Rules');
% paramsToEstimate = {}; 
% initialValues = []; 
% bounds = []; 
% for param_idx = 2:length(Mobj_parameters)
% 
%     add_flag = true; % Whether the parameter should be estimated, default to be true 
%     Mobj_param = Mobj_parameters(param_idx);
% 
%         % Don't estimate the '_F' parameters and just fix them as 1's 
%     if contains(Mobj_param.Name,'_F')
%         add_flag = false;
%     end
% 
%         % Don't estimate AGTPreg_varying 
%     if contains(Mobj_param.Name,'AGTPreg_varying')
%         add_flag = false;
%         % ruleObj = addrule(Model.m1,'AGTPreg_varying=AGTPreg_ON'); % add initial assignment of AGTP regeneration rate
%     end
% 
%     % If parameter is assigned in Model.Rules, no need to estimate 
%     for rule_idx = 1:length(model_rules)
%         rule = model_rules(rule_idx).Rule; 
%         if contains(rule(1:find(rule=='=')),Mobj_param.Name) && ~contains(Mobj_param.Name,'_0')
%             add_flag = false; 
%         end
%     end
% 
%     table_idx = find(strcmp(param_info_table.Name,Mobj_param.Name));
%     LB = table2array(param_info_table(table_idx,'LB')); 
%     UB = table2array(param_info_table(table_idx,'UB')); 
% 
%     %If parameter is fixed (i.e. LB=UB), no need to estimate 
%     if isequal(LB,UB)
%         add_flag = false;
%     end
% 
%     if add_flag
%         paramsToEstimate{end + 1} = strcat('log(',Mobj_param.Name,')');
%         initVal = table2array(param_info_table(table_idx,'InitVal')); 
%         initialValues = [initialValues initVal];
%         bounds = [bounds;[LB,UB]]; 
%     end
% end
% RNA_prod_deg_probObject.Estimated = estimatedInfo(paramsToEstimate,'InitialValue',initialValues,'Bounds',bounds);
% % RNA_deg_probObject.FunctionName = 'lsqcurvefit'; 
% RNA_prod_deg_probObject.Pooled = true; 
% options = optimoptions('lsqnonlin','Display','iter');
% RNA_prod_deg_probObject.Options = options; 
% RNA_prod_deg_probObject.ProgressPlot = 0; 
% 
% % [fitResults,simdataI] = fit(RNA_prod_deg_probObject); % A fit problem object 
% 
% % Run fit 
% all_fitResults = cell(48,1); 
% for iter = 1:48
%     % Modify initial conditions, sample within bounds
%     rng(iter)
%     temp_problemObject = RNA_prod_deg_probObject; 
%     % Use initial parameter values for first run 
%     if ~isequal(iter,1)
%         % Then sample initial parameter values within parameter bounds
%         for est_param_idx = 1:length(temp_problemObject.Estimated)
%             temp_problemObject.Estimated(est_param_idx).InitialValue =  temp_problemObject.Estimated(est_param_idx).Bounds(1) +...
%                 (temp_problemObject.Estimated(est_param_idx).Bounds(2) -  temp_problemObject.Estimated(est_param_idx).Bounds(1)) * rand;
%         end
%     end
%     [fitResults,simdataI] = fit(temp_problemObject); % A fit problem object 
%     all_fitResults{iter,1} = fitResults; 
% end
% 
% all_SSE = nan(48,1);
% for iter = 1:48
%     fitResult = all_fitResults{iter,1}; 
%     all_SSE(iter) = fitResult.SSE;
% end
% save('param_est_run_save/20231019_mRNA_prod_deg_model_fitting.mat');

%% Build a Simbiology model for mRNA degradation & parameter fitting 
%     % Load degradation data. For Broc, the updated data is read on gain =
%     % 70 & gain = 75 
% gain_70_timeVec = readtable('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/processed_data/degradation_curve_gain_70.xlsx','Sheet','timeVec');
% gain_70_RNA_time_course = readtable('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/processed_data/degradation_curve_gain_70.xlsx','Sheet','degradation');
% gain_75_timeVec = readtable('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/processed_data/degradation_curve_gain_75.xlsx','Sheet','timeVec');
% gain_75_RNA_time_course = readtable('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/processed_data/degradation_curve_gain_75.xlsx','Sheet','degradation');
%         % Get initial concentration 
% RNA_init_conc_gain70 = gain_70_RNA_time_course{1,2:9}; 
% RNA_init_conc_gain75 = gain_75_RNA_time_course{1,2:8}; 
%         % Get timeVec and time-course 
% RNA_timeVec_gain70 = gain_70_timeVec{:,2}; 
% RNA_timeVec_gain75 = gain_75_timeVec{:,2}; 
% RNA_fluo_time_course_gain70 = gain_70_RNA_time_course{2:end,11:18}; 
% RNA_fluo_time_course_gain75 = gain_75_RNA_time_course{2:end,10:16}; 
%         % Convert fluorescence to mRNA concentration based on calibration
%         % curve
% % Note: these two should technically give the same mRNA concentrations, but
% % due to experimental noise and fitting erros there's a small difference
% % (~10%) between them. Accounting for uncertainty propogation here would be
% % insane, so let's use gain 70 to fit degradation. 
% RNA_time_course_gain70 = fluo_mRNA_conversion(RNA_fluo_time_course_gain70,[],70); 
% 
%     % Add a time delay to timeVec and data
% RNA_timeVec_gain70 = RNA_timeVec_gain70 + 4 * 60; % 4-min time delay 
% RNA_timeVec_gain70 = [0;RNA_timeVec_gain70]; % Add a true time zero 
% RNA_time_course_gain70 = [RNA_init_conc_gain70;RNA_time_course_gain70]; % Add init conc for the true time zero 
%     % Process data into a Simbiology readable format 
% Time = repmat(RNA_timeVec_gain70,length(RNA_init_conc_gain70),1);
% mRNA_concentration = RNA_time_course_gain70(:); 
%     % Treat negative values as 0 
% mRNA_concentration(mRNA_concentration < 0) = 0; 
% 
% 
% Group = nan(length(mRNA_concentration),1); 
% mRNA_concentration_add = zeros(length(mRNA_concentration),1); 
% for group_idx = 1:length(RNA_init_conc_gain70)
%     Group((group_idx - 1) * length(RNA_timeVec_gain70) + 1:group_idx * length(RNA_timeVec_gain70)) = group_idx; 
%     mRNA_concentration_add((group_idx - 1) * length(RNA_timeVec_gain70) + 1) = RNA_init_conc_gain70(group_idx); 
% end
% RNA_degradation_data_table = array2table([Group,Time,mRNA_concentration,mRNA_concentration_add], 'VariableNames', {'Group','Time', 'mRNA_concentration','mRNA_concentration_add'});
%     % Build a Simbiology model
% Model = sbioloadproject('simbiology_models/RNA_deg.sbproj'); 
% 
% 
%     % Build a fitProblem struct 
% RNA_deg_probObject = fitproblem(); % Initialize a fitProblem object
%         % Data
% RNA_deg_probObject.Data = groupedData(RNA_degradation_data_table); % define a groupedData object 
% RNA_deg_probObject.Data.Properties.IndependentVariableName = 'Time'; 
% RNA_deg_probObject.Data.Properties.GroupVariableName = 'Group'; 
%         % Model 
% RNA_deg_probObject.Model = Model.m1; 
% RNA_deg_probObject.ResponseMap = '[RNA utrbroc--no_protein] = mRNA_concentration';
%         % Doses
% RNA_deg_probObject.Doses = createDoses(RNA_deg_probObject.Data,{'mRNA_concentration_add'}); 
% for dose_idx = 1:size(RNA_deg_probObject.Doses,1)
%     RNA_deg_probObject.Doses(dose_idx,1).TargetName = 'RNA utrbroc--no_protein';
% end
%         % Estimates 
% param_info_table = readtable('parameters_RNAdeg.xlsx');
% Mobj_parameters = get(Model.m1,'Parameters'); 
% paramsToEstimate = {}; 
% initialValues = []; 
% bounds = []; 
% for param_idx = 2:length(Mobj_parameters)
%     Mobj_param = Mobj_parameters(param_idx);
%     table_idx = find(strcmp(param_info_table.Name,Mobj_param.Name));
%     paramsToEstimate{end + 1} = strcat('log(',Mobj_param.Name,')');
%     initVal = table2array(param_info_table(table_idx,'InitVal')); 
%     initialValues = [initialValues initVal];
%     LB = table2array(param_info_table(table_idx,'LB')); 
%     UB = table2array(param_info_table(table_idx,'UB')); 
%     bounds = [bounds;[LB,UB]]; 
% end
% RNA_deg_probObject.Estimated = estimatedInfo(paramsToEstimate,'InitialValue',initialValues,'Bounds',bounds);
% % RNA_deg_probObject.FunctionName = 'lsqcurvefit'; 
% RNA_deg_probObject.Pooled = true; 
% RNA_deg_probObject.ProgressPlot = 0; 
% 
% 
% % [fitResults,simdataI] = fit(RNA_deg_probObject); % A fit problem object 
% 
% % figure;
% % for conc_idx = 1:5%size(RNA_time_course_gain70,2)
% %     plot(RNA_timeVec_gain70,RNA_time_course_gain70(:,conc_idx))
% %     hold on 
% % end
% 
% % Run fit 
% all_fitResults = cell(48,1); 
% for iter = 1:48
%     % Modify initial conditions, sample within bounds
%     rng(iter)
%     temp_problemObject = RNA_deg_probObject; 
%     % Use initial parameter values for first run 
%     if ~isequal(iter,1)
%         % Then sample initial parameter values within parameter bounds
%         for est_param_idx = 1:length(temp_problemObject.Estimated)
%             temp_problemObject.Estimated(est_param_idx).InitialValue =  temp_problemObject.Estimated(est_param_idx).Bounds(1) +...
%                 (temp_problemObject.Estimated(est_param_idx).Bounds(2) -  temp_problemObject.Estimated(est_param_idx).Bounds(1)) * rand;
%         end
%     end
%     [fitResults,simdataI] = fit(temp_problemObject); % A fit problem object 
%     all_fitResults{iter,1} = fitResults; 
% end
% 
% all_SSE = nan(48,1);
% for iter = 1:48
%     fitResult = all_fitResults{iter,1}; 
%     all_SSE(iter) = fitResult.SSE;
% end

%% Random plotting out PE experimental data 
% load param_est_run_save/20231002_param_est_run1014_25.mat
% grouped_data = problemObject.Data;
% [exp_time,exp_data] = process_crosstalk_data_from_source(grouped_data,'grouped_data');
% plot_time_course_crosstalk(exp_time,exp_data)


%% Test toxin mechanism on RNAP
% load param_est_run_save/20231013_param_est_run1302_RNAP_toxin_mech.mat
%     % Plot out all 48 fitting results 
% for iter = 1:length(all_fitResults)
%     fitResult = all_fitResults{iter,1};
%     simData = fitted(fitResult); 
%     [Time,Data] = process_crosstalk_data_from_source(simData,'simulated_data');
%     subplot(7,7,iter)
%     for conc_idx = 1:length(Data)
%         plot(Time{conc_idx},Data{conc_idx},'LineWidth',1.5)
%         hold on 
%     end
%     title(sprintf('iter #%d',iter))
% end
% sgtitle('RNAP toxin mechanism fitting')
% saveas(gcf,'plots/20231016_RNAP_toxin_mechanism_all_fitting','png'); 
%% Confirm that RNA degradation/futile transcription is the cause of decrease 
%     % Experiment #1: Turn off RNA degradation. If hypothesis true, should
%     % not see decrease in expression. 
% load param_est_run_save/20230813_param_est_run0807_PE_T7_strong_no_empty.mat
% fitResult_oi = all_fitResults{36,1}; 
% simulated_data_fit = fitted(fitResult_oi); 
% 
% all_species_name = {problemObject.Model.Species.Name}; 
% sfGFP_RNA_species_idx = contains(all_species_name,'RNA utrGFP--sfGFP') & ~contains(all_species_name,'RNase'); 
% sfGFP_RNA_species_track = all_species_name(sfGFP_RNA_species_idx); 
% additional_track_species = [{'protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA'},sfGFP_RNA_species_track]; 
% resource_names = {'protein sfGFP','protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA','toxin','mRNA sfGFP'}; 
% 
% simFun_ctrl = create_simFun_from_problemObject(problemObject,additional_track_species);
% dosing_information_ctrl = create_dosing_info_from_problemObject(problemObject);
% params = [fitResult_oi.ParameterEstimates.Estimate];
% tStart = 0; 
% tEnd = 21600;
% [ctrl_Time,ctrl_Data] = simFun_ctrl(params',tEnd,dosing_information_ctrl,tStart:tEnd);
% 
%     % Make a copy of the problemObject and revise model in the copy 
% temp_problemObject = problemObject;
% temp_model = temp_problemObject.Model; 
% temp_model_all_rxn = get(temp_model,"Reactions");
% for rxn_idx = 1:length(temp_model_all_rxn)
%     temp_model_rxn = temp_model_all_rxn(rxn_idx); 
%         % If RNase is involved in the reaction, set the reaction to be
%         % inactive 
%     rxn_string = temp_model_rxn.Reaction;
%     if contains(rxn_string,'RNase')
%         set(temp_model_rxn,"Active",false);
%     end
% end
%     % Identify all degradation reactions and set status to inactive 
% simFun = create_simFun_from_problemObject(temp_problemObject,additional_track_species); 
% dosing_information = create_dosing_info_from_problemObject(temp_problemObject); 
% [Time,Data] = simFun(params',tEnd,dosing_information,tStart:tEnd);
% 
% % Plot out resource at 15 nM and 30 nM 
% figure; 
% ctrl_Time_oi = ctrl_Time{6,1}; 
% ctrl_Data_oi = ctrl_Data{6,1}; 
% Time_oi = Time{6,1};
% Data_oi = Data{6,1}; 
% figure;
% for resource_idx = 1:length(resource_names)
%     subplot(4,4,resource_idx)
%     if ~isequal(resource_idx,length(resource_names))
%         plot(ctrl_Time_oi,ctrl_Data_oi(:,resource_idx),'LineWidth',1.5)
%     else
%         plot(ctrl_Time_oi,sum(ctrl_Data_oi(:,resource_idx:end),2),'LineWidth',1.5)
%     end
%     hold on 
%     if ~isequal(resource_idx,length(resource_names))
%         plot(Time_oi,Data_oi(:,resource_idx),'LineWidth',1.5)
%     else
%         plot(Time_oi,sum(Data_oi(:,resource_idx:end),2),'LineWidth',1.5)
%     end
%     title(resource_names{resource_idx})
% 
% end
% legend('with RNA degradation','without RNA degradation')
% 
% figure;
% ctrl_Time_oi = ctrl_Time{7,1}; 
% ctrl_Data_oi = ctrl_Data{7,1}; 
% Time_oi = Time{7,1};
% Data_oi = Data{7,1}; 

%% Test toxin tanh mechanism with RNAP & with AGTP deg time fixed 
% load param_est_run_save/20231011_param_est_run1557_RNAP_toxin_mech.mat
% Plot out predicted for each fitResult
% figure; 
% for i = 1:length(all_fitResults)
%     fitResult = all_fitResults{i,1}; 
%     simData = fitted(fitResult); 
%     [Time,Data] = process_crosstalk_data_from_source(simData,'simulated_data'); 
%     subplot(7,7,i)
    % for data_idx = 1:length(Data)
    %     plot(Time{data_idx},Data{data_idx},'LineWidth',1.5)
    %     hold on 
    % 
    % end
% end
% legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')

% % Get the SSEs 
% all_SSE = nan(length(all_fitResults),1); 
% for i = 1:length(all_fitResults)
%    fitResult = all_fitResults{i,1}; 
%    all_SSE(i) = fitResult.SSE;
% end

% Select the one that captures the decrease in expression/with min SSE
% fitResult_oi = all_fitResults{14,1}; 
% 
% all_species_name = {problemObject.Model.Species.Name}; 
% sfGFP_RNA_species_idx = contains(all_species_name,'RNA utrGFP--sfGFP') & ~contains(all_species_name,'RNase'); 
% sfGFP_RNA_species_track = all_species_name(sfGFP_RNA_species_idx); 
% resource_oi = [{'protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA','toxin'},sfGFP_RNA_species_track]; 
% resource_names = {'protein sfGFP','protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA','toxin','mRNA sfGFP'}; 
% 
% simFunction = create_simFun_from_problemObject(problemObject,resource_oi); 
% dosing_information = create_dosing_info_from_problemObject(problemObject);
% tStart = 0;
% tEnd = 21600; 
% params = [fitResult_oi.ParameterEstimates.Estimate];
% param_names = [fitResult_oi.ParameterEstimates.Name];
% 
% modify_param_name_list = {'AGTP_deg_time_0'}; 
% modify_param_val_list = []; 
% 
% modified_params = modify_params(param_names,params,modify_param_name_list,modify_param_val_list); 
% 
% [Time,Data] = simFunction(modified_params',tEnd,dosing_information,tStart:tEnd);
% 
% figure; 
% for species_idx = 1:length(resource_names)
%     subplot(3,4,species_idx)
%     for data_idx = 1:length(Data)
%         Time_single = Time{data_idx,1};
%         Data_single = Data{data_idx,1}; 
%         if ~isequal(species_idx,length(resource_names))
%             plot(Time_single,Data_single(:,species_idx),'LineWidth',1.5)
%         else
%             plot(Time_single,sum(Data_single(:,species_idx:end)),'LineWidth',1.5)
%         end
%         hold on 
%     end
%     title(resource_names{species_idx})
% end

% group_description = {'PE_T7_strong_no_empty'};
% param_info_path = 'benchmark_model/parameters_test_PE.xlsx';
% 
% weight = []; 
% 
% problemObject_RNAP = setProblemObject_PE_var1_toxRNAP(group_description,param_info_path,weight); 
% problemObject_AGT     Pdegtime = setProblemObject_PE(group_description,param_info_path,weight); 
% 
% wrapper_fit_to_data_w_probObject(problemObject_RNAP,'RNAP_toxin_mech')
% wrapper_fit_to_data_w_probObject(problemObject_AGTPdegtime,'AGTPdegtime_toxin_mech')

%% Determine the order of sequential parameter fitting 
    % run setProblemObject_PE to get both number of parameters to be
    % estimated and number of datapoints 
% param_info_path = '/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/benchmark_model/parameters_test_PE.xlsx';
% 
% round0_group_description  = {'sigma70_strong_no_empty'}; 
% round0_problemObject = setProblemObject_PE(round0_group_description,param_info_path,[]); 
% round0_num_params = length(round0_problemObject.Estimated); 
% round0_num_species = length(round0_problemObject.Model.Species); 
% round0_num_datapoints = numel(round0_problemObject.Data); 

% round1_group_description = {'T7_strong_no_empty','T7_weak_no_empty','sigma70_strong_no_empty'};
% round1_problemObject = setProblemObject_PE(round1_group_description,param_info_path,[]); 
% round1_num_params = length(round1_problemObject.Estimated);
% round1_num_species = length(round1_problemObject.Model.Species); 
% round1_num_datapoints = numel(round1_problemObject.Data); 
% 
% round2_group_description_option1 = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'}; 
% round2_problemObject_option1 = setProblemObject_PE(round2_group_description_option1,param_info_path,[]); 
% round2_num_species_option1 = length(round2_problemObject_option1.Model.Species); 
% round2_num_params_option1 = length(round2_problemObject_option1.Estimated); 
% round2_num_datapoints_option1 = numel(round2_problemObject_option1.Data); 
% round2_param_increase_option1 = round2_num_params_option1 - round1_num_params;
% round2_datapoint_gain_option1 = round2_num_datapoints_option1 - round1_num_datapoints; 
% 
% round2_group_description_option2 = {'T7_strong_no_empty','T7_weak_no_empty','sigma70_strong_no_empty',...
%     'PE_T7_strong_no_empty','PE_T7_weak_no_empty','PE_sigma70_strong_no_empty','PE_sigma70_weak_no_empty'};
% round2_problemObject_option2 = setProblemObject_PE(round2_group_description_option2,param_info_path,[]); 
% round2_num_species_option2 = length(round2_problemObject_option2.Model.Species); 
% round2_num_params_option2 = length(round2_problemObject_option2.Estimated); 
% round2_num_datapoints_option2 = numel(round2_problemObject_option2.Data); 
% round2_param_increase_option2 = round2_num_params_option2 - round1_num_params;
% round2_datapoint_gain_option2 = round2_num_datapoints_option2 - round1_num_datapoints; 
% 
% % How many parameters are added when empty is added?          
% 
% group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70',...
%     'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% 
% problemObject = setProblemObject_PE(group_description,param_info_path,[]); 
% num_params = length(problemObject.Estimated); 
% num_datapoints = numel(problemObject.Data); 

%% Check fitting results for crosstalk ratio fitting for all
% load param_est_run_save/20231004_param_est_run1246.mat
% all_resnorm = nan(length(all_fitResults),1); 
% for i = 1:length(all_fitResults)
%     fitResult = all_fitResults{i,1}; 
%     all_resnorm(i) = fitResult.resnorm; 
% end
% [opt_resnorm,opt_idx] = min(all_resnorm); 
% 
% % Construct problemObject
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = '/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/benchmark_model/parameters_test_PE.xlsx';
% problemObject = setProblemObject_PE(group_description,param_info_path,[]);
% 
% % dosing info
% reporter_prom_name_list = {'DNA pT7--utrGFP--sfGFP','DNA pT773--utrGFP--sfGFP','DNA pJ23119--utrGFP--sfGFP','DNA pJ23105--utrGFP--sfGFP'}; 
% conc_vec = [0.5,1,2.5,5,10,15,30]; 
% empty_conc = 15; 
% num_promotor = length(reporter_prom_name_list); 
% num_conc = 7; 
% tStart = 0; 
% tEnd = 21600; 
% 
% % Create dosing 
% problem_object_dosing_information = create_dosing_info_from_problemObject(problemObject); 
%     % Temp sanity check
% all_species_name = {problemObject.Model.Species.Name}; 
% sfGFP_RNA_species_idx = contains(all_species_name,'RNA utrGFP--sfGFP') & ~contains(all_species_name,'RNase'); 
% sfGFP_RNA_species_track = all_species_name(sfGFP_RNA_species_idx); 
% additional_track_species = [{'protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA'},sfGFP_RNA_species_track]; 
% problemObject_simFunction = create_simFun_from_problemObject(problemObject,additional_track_species); 
% 
% opt_fitResult = all_fitResults{opt_idx,1}; 
% param = opt_fitResult.estimated_params; 
% param_names = [simFunction.Parameters.Name];
% modify_param_name_list = {'TXTL_PT7_RNAPbound_R','TXTL_PT773_RNAPbound_R','TXTL_PJ23119_RNAPbound_R',...
%     'TXTL_PJ23105_RNAPbound_R','TXTL_PkanR_RNAPbound_R','t7RNAP_0','RNAP_0'}; 
% % modify_param_val_list = [0.109336947514019,1.13876756173351,1.14385385734131,...
% %     11.5078292516291,1.12210428662756]; 
% modify_param_val_list = [0.000109336947514019,1.13876756173351,0.114385385734131,...
%     1100000,1.12210428662756,8.7,50]; 
% modified_params = modify_params(param_names,param,modify_param_name_list,modify_param_val_list); 
% 
% % [Time,Data] = simFunction(param',tEnd,dosing_information,tStart:tEnd);
% [problemObject_Time,problemObject_Data] = problemObject_simFunction(modified_params',tEnd,problem_object_dosing_information,tStart:tEnd); 
% 
% experimental_data = problemObject.Data;
% [all_experimental_time,all_experimental_data] = process_crosstalk_data_from_source(experimental_data,'grouped_data'); 
% 
% all_colors = {'g',[0.5,0.5,0.5],'r','b'};
%     figure;
%     for conc_idx = 1:7 % Assume 7 sets of data are always used 
%         plasmid_conc = conc_vec(conc_idx); 
%         for prom_idx = 1:4
%             plot_idx = (conc_idx - 1) * 4 + prom_idx;
%             subplot(7,4,plot_idx)
% 
%             for empty_idx = 1:length(unique_empty_description)
% 
%                 data_idx = (prom_idx - 1) * 28 + (empty_idx - 1) * length(conc_vec) + conc_idx; 
%                 simulated_data = problemObject_Data{data_idx,1}; 
%                 plot(problemObject_Time{data_idx,1},simulated_data(:,1),'LineWidth',1.5,'LineStyle','--','Color',all_colors{empty_idx})
%                 hold on
%                 plot(all_experimental_time{data_idx,1},all_experimental_data{data_idx,1},'LineWidth',1.5,'LineStyle','-','Color',all_colors{empty_idx})
%                 hold on 
%                 if isequal(prom_idx,1)
%                     ylabel(sprintf('%.1f nM',plasmid_conc))
%                 end
%                 if isequal(conc_idx,1)
%                     title(promotors{prom_idx})
%                 end
% 
%             end
%         end
%     end
% 
% line1 = plot(nan,nan,'Color','g','LineStyle','-','DisplayName','no empty experimental');
% line2 = plot(nan,nan,'Color','g','LineStyle','--','DisplayName','no empty simulated');
% line3 = plot(nan,nan,'Color',[0.5,0.5,0.5],'LineStyle','-','DisplayName','empty experimental');
% line4 = plot(nan,nan,'Color',[0.5,0.5,0.5],'LineStyle','--','DisplayName','empty simulated');
% line5 = plot(nan,nan,'Color','r','LineStyle','-','DisplayName','empty T7 experimental');
% line6 = plot(nan,nan,'Color','r','LineStyle','--','DisplayName','empty T7 simulated');
% line7 = plot(nan,nan,'Color','b','LineStyle','-','DisplayName','empty sigma70 experimental');
% line8 = plot(nan,nan,'Color','b','LineStyle','--','DisplayName','empty sigma70 simulated');
% legend([line1,line2,line3,line4,line5,line6,line7,line8])
% % legend([line2,line4,line6,line8])
% 
% % 
% % tStart = 0;
% % tEnd = 21600; 
% % for iter = 1:48
% %     fitResult = all_fitResults{iter,1}; 
% %     param = fitResult.estimated_params; 
% %     [Time,Data] = simFunction(param',tEnd,dosing_information,tStart:tEnd);
% crosstalk_ratio = calculate_crosstalk_ratio_v2(problemObject_Time,problemObject_Data,num_conc,num_promotor,'PE');
% %     if isequal(iter,18) || isequal(iter,47)
% figure; 
% for prom_idx = 1:length(reporter_prom_name_list)
%     subplot(2,2,prom_idx)
%     prom_crosstalk_ratio = crosstalk_ratio(:,prom_idx); 
%     plot(conc_vec,prom_crosstalk_ratio{1},'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%     hold on 
%     plot(conc_vec,prom_crosstalk_ratio{2},'LineWidth',1.5,'Color','r')
%     plot(conc_vec,prom_crosstalk_ratio{3},'LineWidth',1.5,'Color','b')
%     title(promotors{prom_idx})
% end
% end
%     % figure; 
%     % for conc_idx = 1:length(conc_vec)
%     %     subplot(3,3,conc_idx)
%     %     plot(Time{conc_idx,1},Data{conc_idx,1},'LineWidth',1.5,'Color','g')
%     %     hold on 
%     %     plot(Time{conc_idx + 1 * length(conc_vec),1},Data{conc_idx + 1 * length(conc_vec),1},'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%     %     plot(Time{conc_idx + 2 * length(conc_vec),1},Data{conc_idx + 2 * length(conc_vec),1},'LineWidth',1.5,'Color','r')
%     %     plot(Time{conc_idx + 3 * length(conc_vec),1},Data{conc_idx + 3 * length(conc_vec),1},'LineWidth',1.5,'Color','b')
%     %     title(sprintf('%.1f nM',conc_vec(conc_idx)))
%     % end
%     % sgtitle(sprintf('PE all crosstalk - ratio as obj fun (iter #%d)',iter))
% end

%% Check processed large-scale optimization run results 
% timestamp_confusion = {'1531','1532','1533','1534'}; 
% result_fileName_stem = 'param_est_run_save/20231102_param_est_run'; 
% num_conc = 7;
% num_promotor = 4; 
% sz = [4800 5];
% varTypes = ["string","double","double","cell","cell"];
% varNames = ["Label","SSE","LogLikelihood","CrosstalkRatios","ParameterEstimates"];
% optimization_run_table = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 
% for run_idx = 1:100
%     result_fileName_suffix = sprintf('_%d.mat',run_idx); 
%     for timeStamp_idx = 1:length(timestamp_confusion)
%         timestamp = timestamp_confusion{timeStamp_idx}; 
%         result_fileName = strcat(result_fileName_stem,timestamp,result_fileName_suffix); 
%         try
%             result_file = load(result_fileName); 
%             break
%         catch
%             continue
%         end
%     end
%     try 
%     all_fitResults = result_file.all_fitResults; 
%     problemObject = result_file.problemObject; 
%     % For each fitted result, calculate both deviation in absolute values
%     % and in crosstalk ratios (i.e. whether both positive & negative
%     % crosstalk captured) 
%     for iter = 1:length(all_fitResults)
%         fitResult = all_fitResults{iter};
%         all_simData = fitted(fitResult); 
%         [Time,Data] = process_crosstalk_data_from_source(all_simData,'simulated_data'); 
%         crosstalk_ratio = calculate_crosstalk_ratio_v2(Time,Data,num_conc,num_promotor,'PE'); 
%         label = sprintf('run%d_iter%d',run_idx,iter); 
%         SSE = fitResult.SSE; 
%         LogLikelihood = fitResult.LogLikelihood; 
%         ParameterEstimates = [fitResult.ParameterEstimates.Estimate]; 
% 
%         optimization_run_table((run_idx - 1) * length(all_fitResults) + iter,:) = {label,SSE,LogLikelihood,crosstalk_ratio,ParameterEstimates}; 
% 
%     end
%     fprintf('%d',run_idx); 
%     end
% end
% save('param_est_run_save/20231106_large_scale_run.mat','optimization_run_table')
%% Analyze large-scale parameter fitting results 
    % Initial processing -> get metrics and best fit 
% timestamp_confusion = {'1531','1532','1533','1534'}; 
% result_fileName_stem = 'param_est_run_save/20231102_param_est_run'; 
% for iter = 1:100
%     result_fileName_suffix = sprintf('_%d.mat',iter); 
%     for timeStamp_idx = 1:length(timestamp_confusion)
%         timestamp = timestamp_confusion{timeStamp_idx}; 
%         result_fileName = strcat(result_fileName_stem,timestamp,result_fileName_suffix); 
%         try
%             result_file = load(result_fileName); 
%             break
%         catch
%             continue
%         end
%     end
%     all_fitResults = result_file.all_fitResults; 
%     problemObject = result_file.problemObject; 
%     fit_result_path = result_fileName(1:end-4); 
%     metric_name = 'SSE'; 
%     plot_fitting = false;
%     plot_param_dist = false; 
%     wrapper_analyze_fit_result_v2(all_fitResults,problemObject,fit_result_path,metric_name,plot_fitting,plot_param_dist) 
% 
% 
% 
% end

%% Fit crosstalk ratio for sigma70 weak 
% group_description = {'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = '/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/benchmark_model/parameters_test_PE.xlsx';
% time = datetime; 
% save_path_id = ''; 
% save_path_stem = sprintf('param_est_run_save/2023%02d%02d_param_est_run%02d%02d_%s',time.Month,time.Day,time.Hour,time.Minute,save_path_id); 
% wrapper_fit_crosstalk_ratio(group_description,param_info_path,save_path_stem)


%% Check crosstalk in weak promotors 

% fit_result_path = '20230813_param_est_run0807_sigma70_weak_PE'; 
% metric_name = 'SSE'; 
% plot_fitting = true; 
% plot_param_dist = true; 
% wrapper_analyze_fit_result(fit_result_path,metric_name,plot_fitting,plot_param_dist) 

% load param_est_run_save/20230813_param_est_run0807_sigma70_weak_PE_analysis.mat
% 
% % Calculate crosstalk ratio for each fitResult
% % ,additional_track_species
% simFunction = create_simFun_from_problemObject(problemObject); 
% tStart = 0;
% tEnd = 21600; 
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% all_crosstalk_ratio = []; 
% for iter = 1:length(metric_summary)
%     eval(sprintf('params = [estimated_param_summary.Estimate_%d];',iter))
%     [Time,Data] = simFunction(params',tEnd,dosing_information,tStart:tEnd);
% 
%     crosstalk_ratio = calculate_crosstalk_ratio_v2(Time,Data,num_conc,1,'PE'); 
%     all_crosstalk_ratio= [all_crosstalk_ratio,crosstalk_ratio]; 
% end
% 
% figure; 
% for iter = 1:length(metric_summary)
%     subplot(7,7,iter)
%     plot(conc_vec,all_crosstalk_ratio{1,iter},'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%     hold on 
%     plot(conc_vec,all_crosstalk_ratio{2,iter},'LineWidth',1.5,'Color','r')
%     plot(conc_vec,all_crosstalk_ratio{3,iter},'LineWidth',1.5,'Color','b')
% end

% % Focus on the opt_fitesul
% % opt_fitResult_file = load('param_est_run_save/20230813_param_est_run0807_sigma70_weak_PE_fitResult46.mat'); 
% opt_fitResult = opt_fitResult_file.fitResults; 
% % simData = fitted(opt_fitResult); 
% all_species_name = {problemObject.Model.Species.Name}; 
% sfGFP_RNA_species_idx = contains(all_species_name,'RNA utrGFP--sfGFP') & ~contains(all_species_name,'RNase'); 
% sfGFP_RNA_species_track = all_species_name(sfGFP_RNA_species_idx); 
% resource_oi = [{'protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA'},sfGFP_RNA_species_track]; 
% simFunction = create_simFun_from_problemObject(problemObject,resource_oi); 
% % Modify parameters 
% params = [opt_fitResult.ParameterEstimates.Estimate];
% param_names = [opt_fitResult.ParameterEstimates.Name];
% modify_param_name_list = {'AGTPdeg_time','RNase_0','TX_elong_glob','TXTL_PJ23105_RNAPbound_R'}; 
% modify_param_val_list = [4000,20,7,1e+05]; 
% modified_params = modify_params(param_names,params,modify_param_name_list,modify_param_val_list); 
% [Time,Data] = simFunction(modified_params',tEnd,dosing_information,tStart:tEnd);
%     % Visualize time-course GFP
% all_color_list = {'g',[0.5,0.5,0.5],'r','b'}; 
% 
% resource_oi_name = {'protein sfGFP*','protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA','RNA sfGFP'}; 
% figure; 
% for comb_idx = 1:length(unique_empty_description)
%     for conc_idx = 1:num_conc
% 
%         conc = conc_vec(conc_idx); 
%         simData_idx = (comb_idx - 1) * num_conc + conc_idx; 
%         timeVec = Time{simData_idx,1}; 
%         simData_single = Data{simData_idx,1};
% 
%         for resource_idx = 1:length(resource_oi_name)
%             subplot(7,length(resource_oi_name),(conc_idx - 1) * length(resource_oi_name) + resource_idx)
%             if isequal(resource_idx,length(resource_oi_name))
%                 plot(timeVec,sum(simData_single(:,resource_idx:end),2),'LineWidth',1.5,'Color',all_color_list{comb_idx})
%             else
%                 plot(timeVec,simData_single(:,resource_idx),'LineWidth',1.5,'Color',all_color_list{comb_idx})
%             end
%             hold on 
%             if isequal(conc_idx,1)
%                 title(resource_oi_name{resource_idx});
%             end
%             if isequal(resource_idx,1)
%                 ylabel(sprintf('%.1f nM',conc_vec(conc_idx)))
%             end
%         end
%     end
% end




%% Check 20230921_0948 
% load param_est_run_save/20230921_param_est_run0948_.mat
% % observable_name = 'protein sfGFP*'; 
% metric_summary = nan(48,1); 
% for iter = 1:48
%     fitResult = all_fitResults{iter,1}; 
%     if ~isempty(fitResult)
%         metric_summary(iter) = fitResult.SSE; 
%     end
% end
% [~,opt_idx] = min(metric_summary); 
% opt_fitResult = all_fitResults{opt_idx,1}; 

% opt_simData = fitted(opt_fitResult); 
% figure; 
% for conc_idx = 1:length(conc_vec)
%     timeVec = Time{conc_idx,1}; 
%     dataVec = Data{conc_idx,1}; 
%     % plot(Time{conc_idx,1},Data{conc_idx,1},'LineWidth',1.5)
%     plot(timeVec,dataVec(:,1),'LineWidth',1.5)
%     hold on
% end
% legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')
% all_species_name = {problemObject.Model.Species.Name}; 
% sfGFP_RNA_species_idx = contains(all_species_name,'RNA utrGFP--sfGFP') & ~contains(all_species_name,'RNase'); 
% sfGFP_RNA_species_track = all_species_name(sfGFP_RNA_species_idx); 
% additional_track_species = [{'AGTP','AA','protein kanR','RNAP','Ribo'},sfGFP_RNA_species_track]; 
% simFunction = create_simFun_from_problemObject(problemObject,additional_track_species); 
% tStart = 0;
% tEnd = 21600; 
% params = [opt_fitResult.ParameterEstimates.Estimate];
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
% [Time,Data] = simFunction(params',tEnd,dosing_information,tStart:tEnd);

%% Check fitResult #12 & #36 for T7_strong_no_empty fit 
    % For the decrease in PE at high plasmid concentration 
% load param_est_run_save/20230813_param_est_run0807_PE_T7_strong_no_empty.mat
% fitResult_oi = all_fitResults{36,1};
% fitted_data_oi = fitted(fitResult_oi);

% Get all species relevant to sfGFP transcription 
% all_species_name = {problemObject.Model.Species.Name}; 
% relevant_species_names = all_species_name(contains(all_species_name,'sfGFP')); 
% [timeVec,sfGFP_relevant_tc] = get_species_time_course(fitted_data_oi(7,1),relevant_species_names); 

% % Filter out time-course that is always close to zero 
% keep_idx = ~all(sfGFP_relevant_tc < 1e-08,1);
% relevant_species_names_nonzero = relevant_species_names(keep_idx);
% sfGFP_relevant_tc_nonzero = sfGFP_relevant_tc(:,keep_idx); 
% 
% 
% % Plot all relevant sfGFP species 
% figure;
% for subplot_idx = 1:size(sfGFP_relevant_tc_nonzero,2)
%     subplot(ceil(sqrt(size(sfGFP_relevant_tc_nonzero,2))),ceil(sqrt(size(sfGFP_relevant_tc_nonzero,2))),subplot_idx)
%     plot(timeVec,sfGFP_relevant_tc_nonzero(:,subplot_idx),'LineWidth',1.5)
%     hold on 
%     title(relevant_species_names_nonzero{subplot_idx})
% end

% Plot CUTP
% resource_name_list = {'AGTP','CUTP','AA','RNAP','t7RNAP','Ribo'};
% for conc_idx = 1:length(conc_vec)
% 
%     [timeVec,resource_tc] = get_species_time_course(fitted_data_oi(conc_idx,1),resource_name_list); 
%     for subplot_idx = 1:size(resource_tc,2)
%         subplot(2,3,subplot_idx)
%         plot(timeVec,resource_tc(:,subplot_idx),'LineWidth',1.5)
%         hold on 
%         title(resource_name_list{subplot_idx}) 
%     end
% 
% end
% Play with parameters in # 36
%     % Create sim Function 
% all_species_name = {problemObject.Model.Species.Name}; 
% sfGFP_RNA_species_idx = contains(all_species_name,'RNA utrGFP--sfGFP') & ~contains(all_species_name,'RNase'); 
% sfGFP_RNA_species_track = all_species_name(sfGFP_RNA_species_idx); 
% additional_track_species = [{'AGTP','CUTP','AA','protein kanR','RNAP','Ribo'},sfGFP_RNA_species_track]; 
% simFunction = create_simFun_from_problemObject(problemObject,additional_track_species); 
% tStart = 0;
% tEnd = 21600; 
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
%     % Increasing the dose to see the impact of CUTP resource competition
% modified_dosing_information = dosing_information;
% modified_dosing_information{1,4}.Amount = 20; 
% modified_dosing_information{1,7}.Amount = 20; 
% modified_dosing_information{2,4}.Amount = 25; 
% modified_dosing_information{2,7}.Amount = 25; 
% modified_dosing_information{3,4}.Amount = 100; 
% modified_dosing_information{3,7}.Amount = 100; 
% modified_dosing_information{4,4}.Amount = 50; 
% modified_dosing_information{4,7}.Amount = 50;
% 
% 
% %     % Change parametera
% params = [fitResult_oi.ParameterEstimates.Estimate];
% %     param_names = [fitResult_oi.ParameterEstimates.Name];
% %     modify_param_name_list = {'AGTPdeg_time'}; 
% %     modify_param_val_list = 3900;
% %     modified_params = modify_params(param_names,params,modify_param_name_list,modify_param_val_list);
% [Time,Data] = simFunction(params',tEnd,modified_dosing_information,tStart:tEnd);
% % 
% %     % Check AGTP consumption situation here 
% %     % Let's directly compare all resource utilization for 0.5nM, 5nM, 10nM &
% % % 30nM 
% conc_vec_oi = [20,25,40,50,10,15,30]; 
% figure; 
% for conc_idx = 1:length(conc_vec_oi)
%     timeVec = Time{conc_idx,1}; 
%     dataVec = Data{conc_idx,1}; 
%     % plot(Time{conc_idx,1},Data{conc_idx,1},'LineWidth',1.5)
%     plot(timeVec,dataVec(:,1),'LineWidth',1.5)
%     hold on
% end
% legend('20nM','25nM','40nM','50nM','10nM','15nM','30nM')
% % legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')
% % 
% figure; 
% for conc_idx = 1:length(conc_vec_oi)
%      % simData = fitted_data_oi(conc_idx,1); 
%      timeVec = Time{conc_idx,1}; 
%      dataVec = Data{conc_idx,1}; 
% 
%      % sfGFP_idx = strcmp(simData.DataNames,'protein sfGFP*');
%      % sfGFP_tc = simData.Data(:,sfGFP_idx); 
%      % plot(simData.Time,sfGFP_tc,'LineWidth',1.5)
%      % hold on 
% 
%      % AGTP_idx = strcmp(simData.DataNames,'AGTP'); 
%      % AA_idx = strcmp(simData.DataNames,'AA'); 
%      % protein_kanR_idx = strcmp(simData.DataNames,'protein kanR');
%      % RNAP_idx = strcmp(simData.DataNames,'RNAP'); 
%      % ribosome_idx = strcmp(simData.DataNames,'Ribo'); 
%      % sfGFP_mRNA_idx = contains(simData.DataNames,'RNA utrGFP--sfGFP') & ~contains(simData.DataNames,'RNase'); 
% 
%      % AGTP_time_course = simData.Data(:,AGTP_idx); 
%      % AA_time_course = simData.Data(:,AA_idx);
%      % protein_kanR_tc = simData.Data(:,protein_kanR_idx); 
%      % RNAP_time_course = simData.Data(:,RNAP_idx); 
%      % Ribo_time_course = simData.Data(:,ribosome_idx); 
%      % sfGFP_mRNA_tc = sum(simData.Data(:,sfGFP_mRNA_idx),2); 
%      AGTP_time_course = dataVec(:,2); 
%      CUTP_time_course = dataVec(:,3); 
%      AA_time_course = dataVec(:,4);
%      protein_kanR_tc = dataVec(:,5); 
%      RNAP_time_course = dataVec(:,6); 
%      Ribo_time_course = dataVec(:,7); 
%      sfGFP_mRNA_tc = sum(dataVec(:,8:end),2); 
% 
%      subplot(3,3,1)
%      plot(timeVec,AGTP_time_course,'LineWidth',1.5)
%      hold on 
%      title('AGTP')
%      subplot(3,3,2)
%      plot(timeVec,CUTP_time_course,'LineWidth',1.5)
%      hold on 
%      title('CUTP')
%      subplot(3,3,3)
%      plot(timeVec,AA_time_course,'LineWidth',1.5)
%      hold on 
%      title('AA')
%      subplot(3,3,4)
%      plot(timeVec,protein_kanR_tc,'LineWidth',1.5)
%      hold on 
%      title('protein kanR')
%      subplot(3,3,5)
%      plot(timeVec,RNAP_time_course,'LineWidth',1.5)
%      hold on 
%      title('RNAP')
%      subplot(3,3,6)
%      plot(timeVec,Ribo_time_course,'LineWidth',1.5)
%      hold on 
%      title('Ribo')
%      subplot(3,3,7)
%      plot(timeVec,sfGFP_mRNA_tc,'LineWidth',1.5)
%      hold on 
%      title('sfGFP mRNA')
% end
% legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')

% Let's directly compare all resource utilization for 0.5nM, 5nM, 10nM &
% 30nM 

%% Take a fitted T7_strong_no_empty result and change AGTP_deg_time
% load param_est_run_save/20230813_param_est_run0807_PE_T7_strong_no_empty.mat
% for iter = 1:length(all_fitResults)
%     fitResult = all_fitResults{iter,1}; 
%     params = [fitResult.ParameterEstimates.Estimate];
% 
%     % Create sim Function 
%     simFunction = create_simFun_from_problemObject(problemObject); 
%     tStart = 0;
%     tEnd = 21600; 
%     dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
%     % Run with originally fitted parameters 
%     [ctrl_Time,ctrl_Data] = simFunction(params',tEnd,dosing_information,tStart:tEnd);
% 
%     % Change AGTP_deg_time to a smaller number and run sim again  
%     param_names = [fitResult.ParameterEstimates.Name];
%     modify_param_name_list = {'AGTPdeg_time'}; 
%     original_param_val = fitResult.ParameterEstimates.Estimate(strcmp(param_names,'AGTPdeg_time')); 
%     modify_param_val_list = original_param_val * 0.5;
%     modified_params = modify_params(param_names,params,modify_param_name_list,modify_param_val_list);
%     [Time,Data] = simFunction(modified_params',tEnd,dosing_information,tStart:tEnd);
% 
%     % Visualize and compare the simulations 
%     figure(iter);
%         % Put control on the left 
%     subplot(1,2,1)
%     for ctrl_data_idx = 1:length(ctrl_Data)
%         plot(ctrl_Time{ctrl_data_idx,1},ctrl_Data{ctrl_data_idx},'LineWidth',1.5)
%         hold on 
%     end
%         % And the modified on the right
%     subplot(1,2,2)
%     for data_idx = 1:length(Data)
%         plot(Time{data_idx,1},Data{data_idx},'LineWidth',1.5)
%         hold on 
%     end
% end

%% Check 0919_1229 rate mechanism results 
% load param_est_run_save/20230919_param_est_run1229_.mat
% metric_summary = nan(length(all_fitResults),1); 
% for iter = 1:length(all_fitResults)
%     fitResult = all_fitResults{iter,1}; 
%     if ~isempty(fitResult)
%         metric_summary(iter) = fitResult.SSE; 
%     end
% 
% end
% [min_SSE,min_idx] = min(metric_summary);
% opt_fitResult = all_fitResults{min_idx,1}; 
% fitted_data = fitted(opt_fitResult); 
% opt_params = [opt_fitResult.ParameterEstimates.Estimate];
% [Time,Data] = process_crosstalk_data_from_source(fitted_data,'simulated_data'); 
% additional_track_species = {'AGTP','toxin'};
% simFunction = create_simFun_from_problemObject(problemObject,additional_track_species); 
% tStart = 0;
% tEnd = 21600; 
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
% 
% [Time,Data] = simFunction(opt_params',tEnd,dosing_information,tStart:tEnd);
% figure; 
% for data_idx = 1:length(fitted_data)
% 
%     timeVec = Time{data_idx,1};
%     dataVec = Data{data_idx,1};
%     plot(timeVec,dataVec(:,2),'LineWidth',1.5)
%     hold on 
% end

%% Sanity check crosstalk ratio as objective function
% crosstalk_ratio_fit_result = load('param_est_run_save/20230831_param_est_run1323.mat'); 
% PE_crosstalk_fit_result = load('param_est_run_save/20230802_param_est_run1352_.mat'); 
% 
% % Get the optimal result for ratio and the worst result for PE_crosstalk
% % fit 
% PE_crosstalk_opt_ratio_fitResult = PE_crosstalk_fit_result.all_fitResults{37,1}; 
% crosstalk_ratio_opt_fitResult = crosstalk_ratio_fit_result.all_fitResults{45,1}; 
% 
% % Calculate crosstalk ratio based on the data
% 
%     % PE_crosstalk
% PE_crosstalk_predicted_data = fitted(PE_crosstalk_opt_ratio_fitResult); 
% [PE_crosstalk_predicted_Time,PE_crosstalk_predicted_Data] = ...
%     process_crosstalk_data_from_source(PE_crosstalk_predicted_data,'simulated_data'); 
% mode = 'PE'; 
% PE_crosstalk_ratio = calculate_crosstalk_ratio_v2(PE_crosstalk_predicted_Time,PE_crosstalk_predicted_Data,num_conc,num_promotor,mode);
% 
%     % Crosstalk ratio
% tStart = 0; 
% tEnd = 21600;
% opt_estimated_params = crosstalk_ratio_opt_fitResult.estimated_params;
% simFunction = crosstalk_ratio_fit_result.simFunction; 
% dosing_information = create_dosing_info_from_problemObject(PE_crosstalk_fit_result.problemObject);
% filtered_dosing_information = dosing_information(:,[7:10,4:6]); 
% [crosstalk_ratio_Time,crosstalk_ratio_Data] = simFunction(opt_estimated_params',tEnd,filtered_dosing_information,tStart:tEnd);
% crosstalk_ratio_obj = calculate_crosstalk_ratio_v2(crosstalk_ratio_Time,crosstalk_ratio_Data,num_conc,num_promotor,mode);
% 
% % Compare the crosstalk ratio
%     % load the actual crosstalk ratio 
% grouped_data = get_data_table(group_description); 
% [crosstalk_data_time,crosstalk_data_data] = ...
%     process_crosstalk_data_from_source(grouped_data,'grouped_data'); 
% exp_crosstalk_ratio = calculate_crosstalk_ratio_v2(crosstalk_data_time,crosstalk_data_data,num_conc,num_promotor,mode);
%     % Calculate residual norm of the two result structs
% PE_fit_crosstalk_ratio_expanded = [PE_crosstalk_ratio{:}];
% crosstalk_ratio_obj_expanded = [crosstalk_ratio_obj{:}];
% exp_crosstalk_ratio_expanded = [exp_crosstalk_ratio{:}];
% PE_fit_error = sum(sum((PE_fit_crosstalk_ratio_expanded - exp_crosstalk_ratio_expanded).^2));
% crosstalk_ratio_obj_error = sum(sum((crosstalk_ratio_obj_expanded - exp_crosstalk_ratio_expanded).^2)); 
% 
% % Use the parameters from PE_crosstalk as the initial point for crosstalk
% % ratio fitting 
% PE_fit_crosstalk_opt_params = [PE_crosstalk_opt_ratio_fitResult.ParameterEstimates.Estimate];
% PE_fit_crosstalk_opt_params = [PE_fit_crosstalk_opt_params(1:10);PE_fit_crosstalk_opt_params(12:end)];
% [sanity_check_time,sanity_check_data] = simFunction(PE_fit_crosstalk_opt_params',tEnd,filtered_dosing_information,tStart:tEnd);
% sanity_check_ratio = calculate_crosstalk_ratio_v2(sanity_check_time,sanity_check_data,num_conc,num_promotor,mode);

%% Check the fitting of crosslk ratio as objective function 
% load param_est_run_save/20230831_param_est_run1323.mat
% all_resnorm = nan(length(crosstalk_ratio_fit_result.all_fitResults),1);
% all_fitted_params = nan(length(crosstalk_ratio_fit_result.all_fitResults{1,1}.estimated_params),length(crosstalk_ratio_fit_result.all_fitResults)); 
% for i = 1:length(crosstalk_ratio_fit_result.all_fitResults)
%     all_resnorm(i) = crosstalk_ratio_fit_result.all_fitResults{i,1}.resnorm;
%     all_fitted_params(:,i) = crosstalk_ratio_fit_result.all_fitResults{i,1}.estimated_params;
% end
% [~,opt_idx] = min(all_resnorm);
% opt_fitResult = crosstalk_ratio_fit_result.all_fitResults{opt_idx,1}; 
% opt_params = all_fitted_params(:,opt_idx); 
% 
% % Compile dosing information 
% reporter_prom_name_list = {'DNA pT7--utrGFP--sfGFP','DNA pT773--utrGFP--sfGFP','DNA pJ23119--utrGFP--sfGFP',...
%         'DNA pJ23105--utrGFP--sfGFP'};
% num_promotor = length(reporter_prom_name_list); 
% num_conc = length(conc_vec); 
% dosing_information_mat = zeros(num_promotor * num_conc * 4,length(simFunction.Dosed.TargetName));
% for reporter_prom_idx = 1:length(reporter_prom_name_list)
%     reporter_prom_name = reporter_prom_name_list{reporter_prom_idx}; 
%     start_idx = (reporter_prom_idx - 1) * length(conc_vec) * 4; 
% 
%     for conc_idx = 1:length(conc_vec)
%         reporter_conc = conc_vec(conc_idx);
%         reporter_idx = strcmp(reporter_prom_name,simFunction.Dosed.TargetName); 
% 
%         dosing_information_mat(start_idx + conc_idx,reporter_idx) = reporter_conc; 
%         dosing_information_mat(start_idx + 1 * length(conc_vec) + conc_idx,reporter_idx) = reporter_conc; 
%         dosing_information_mat(start_idx + 2 * length(conc_vec) + conc_idx,reporter_idx) = reporter_conc; 
%         dosing_information_mat(start_idx + 3 * length(conc_vec) + conc_idx,reporter_idx) = reporter_conc; 
% 
%         kanR_idx = strcmp('DNA pkanR--utrkanR--kanR',simFunction.Dosed.TargetName);
%         kanR_conc = empty_conc + reporter_conc;  
% 
%         dosing_information_mat(start_idx + conc_idx,kanR_idx) = reporter_conc; 
%         dosing_information_mat(start_idx + 1 * length(conc_vec) + conc_idx,kanR_idx) = kanR_conc; 
%         dosing_information_mat(start_idx + 2 * length(conc_vec) + conc_idx,kanR_idx) = kanR_conc; 
%         dosing_information_mat(start_idx + 3 * length(conc_vec) + conc_idx,kanR_idx) = kanR_conc; 
%     end
%     empty_T7_idx = strcmp('DNA pT7--utrempty--no_protein',simFunction.Dosed.TargetName);
%     empty_sigma70_idx = strcmp('DNA pJ23119--utrempty--no_protein',simFunction.Dosed.TargetName);
%     dosing_information_mat(start_idx + 2 * length(conc_vec) + 1:start_idx + 2 * length(conc_vec) + length(conc_vec),empty_T7_idx) = empty_conc; 
%     dosing_information_mat(start_idx + 3 * length(conc_vec) + 1:start_idx + 3 * length(conc_vec) + length(conc_vec),empty_sigma70_idx) = empty_conc; 
% end
% dosing_information = cell(size(dosing_information_mat)); 
% for row_idx = 1:size(dosing_information_mat,1)
%     for col_idx = 1:size(dosing_information_mat,2)
%         dosing_table = array2table([0,dosing_information_mat(row_idx,col_idx)]); 
%         dosing_table.Properties.VariableNames = {'Time','Amount'};
%         dosing_information{row_idx,col_idx} = dosing_table; 
%     end
% end
% 
% tStart = 0; 
% tEnd = 21600; 
% 
% [Time,Data] = simFunction(opt_params',tEnd,dosing_information,tStart:tEnd);
% % Calculate crosstalk ratios from simulated data 
% simulated_crosstalk_ratios = calculate_crosstalk_ratio_simFunction(Time,Data,num_conc,num_promotor,'PE'); 
% 
% % Actual simulated data
% figure;
% count_idx = 1; 
% for conc_idx = 1:length(conc_vec) % Assume 7 sets of data are always used 
%     plasmid_conc = conc_vec(conc_idx); 
%     for prom_idx = 1:length(promotors)
% 
%         plot_idx = (conc_idx - 1) * length(promotors) + prom_idx;
%         subplot(7,length(promotors),plot_idx)
% 
%         promotor_description = promotors{prom_idx};
% 
%         if isequal(conc_idx,1)
%             title(strrep(promotor_description,'_',' '))
%             hold on 
%         end
%         if isequal(conc_idx,7)
%             xlabel('time(s)')
%             hold on 
%         end
%         if isequal(prom_idx,1)
%             ylabel(sprintf('%.1f nM',plasmid_conc))
%             hold on 
%         end
% 
%         for empty_idx = 1:length(unique_empty_description)
%             empty_description = unique_empty_description{empty_idx}; 
% 
%             data_idx = (prom_idx - 1) * length(conc_vec) * 4 + (empty_idx - 1) * length(conc_vec) + conc_idx; 
% 
%             plot(Time{data_idx,1},Data{data_idx,1},'LineWidth',1.5,'LineStyle','--','Color',all_colors{empty_idx})
%             hold on
%             plot(Time{data_idx,1},Data{data_idx,1},'LineWidth',1.5,'LineStyle','-','Color',all_colors{empty_idx})
%             hold on 
% 
%         end
%     end
% end
% 
% % line1 = plot(nan,nan,'Color','g','LineStyle','-','DisplayName','no empty experimental');
% line2 = plot(nan,nan,'Color','g','LineStyle','--','DisplayName','no empty simulated');
% % line3 = plot(nan,nan,'Color',[0.5,0.5,0.5],'LineStyle','-','DisplayName','empty experimental');
% line4 = plot(nan,nan,'Color',[0.5,0.5,0.5],'LineStyle','--','DisplayName','empty simulated');
% % line5 = plot(nan,nan,'Color','r','LineStyle','-','DisplayName','empty T7 experimental');
% line6 = plot(nan,nan,'Color','r','LineStyle','--','DisplayName','empty T7 simulated');
% % line7 = plot(nan,nan,'Color','b','LineStyle','-','DisplayName','empty sigma70 experimental');
% line8 = plot(nan,nan,'Color','b','LineStyle','--','DisplayName','empty sigma70 simulated');
% % legend([line1,line2,line3,line4,line5,line6,line7,line8])
% legend([line2,line4,line6,line8])


% % Crosstalk ratio predicted by opt_fitResult
% figure; 
% color_list = {[0.5,0.5,0.5],'r','b'};
% for prom_idx = 1:size(simulated_crosstalk_ratios,2)
%     subplot(2,2,prom_idx)
%     for comp_idx = 1:size(simulated_crosstalk_ratios,1)
%         crosstalk_ratio = simulated_crosstalk_ratios{comp_idx,prom_idx};
%         plot(conc_vec,crosstalk_ratio,'Color',color_list{comp_idx},'LineWidth',0.5)
%         hold on 
%         title(promotors{prom_idx})
%     end
% end
% sgtitle('Crosstalk ratio fitting results')
% 
% % Distribution of residues 
% figure;
% bar(all_resnorm)

%% Test process different data type for consolidating calculate_crosstalk_ratio 
    % grouped_data
%     group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%         'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%         'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%         'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
%     % grouped_data = get_data_table(group_description);
%     % source = 'grouped_data';
%     % % simulated_data
%     % source = 'simulated_data';
%     % % sim_function
%     % source = 'sim_function'; 
% [Time,Data] = process_crosstalk_data_from_source(grouped_data,source);

% function grouped_data = get_data_table(group_description)
%     %% Load data 
%         % load both transcription and protein expression data 
%     tx_data_file = load('data_structures/simbio_data_table_updated_FP.mat'); 
%     PE_data_file = load('data_structures/simbio_data_table_PE_updated_FP.mat');
%     tx_data_table = tx_data_file.all_data_table; 
%     PE_data_table = PE_data_file.all_data_table; 
%     tx_data_description = tx_data_file.variable_name_stem_list; 
%     PE_data_description = PE_data_file.variable_name_stem_list; 
% 
%         % Select relevant data 
%     [tx_group_number,tx_selected_table] = get_relevant_data(tx_data_table,tx_data_description,group_description); 
%     [PE_group_number,PE_selected_table] = get_relevant_data(PE_data_table,PE_data_description,group_description); 
% 
%         % Convert fluorescence to concentration and add concentration to
%         % data table 
%         %%%%%% mRNA %%%%%%
%     if ~isempty(tx_selected_table)
%         split_point_idx_list = find(tx_selected_table.gain==tx_selected_table.gain(1)); 
%         split_point_idx = split_point_idx_list(end); 
%         if isequal(length(split_point_idx_list),height(tx_selected_table))% if there's only one gain 
%             mRNA_concentration = [fluo_mRNA_conversion(tx_selected_table.fluorescence,[],tx_selected_table.gain(1))];
%         else
%                 % Applies to when there are 2 gains 
%             mRNA_concentration = [fluo_mRNA_conversion(tx_selected_table.fluorescence(1:split_point_idx),[],tx_selected_table.gain(1));...
%                 fluo_mRNA_conversion(tx_selected_table.fluorescence(split_point_idx+1:end),[],tx_selected_table.gain(split_point_idx + 1))];
%         end
%         tx_selected_table = addvars(tx_selected_table,mRNA_concentration,'NewVariableNames','mRNA_concentration'); 
%             % Replace negative values with nan 
%         tx_selected_table = remove_neg_val(tx_selected_table);
%             % Find new split after negative value removal
%         new_split_point_idx_list = find(tx_selected_table.gain==tx_selected_table.gain(1)); 
%         new_split_point_idx = new_split_point_idx_list(end); 
%         if isequal(length(new_split_point_idx_list),height(tx_selected_table))
%             mRNA_CI_lb = [fluo_mRNA_conversion(tx_selected_table.CI_lb,[],tx_selected_table.gain(1))];
%             mRNA_CI_ub = [fluo_mRNA_conversion(tx_selected_table.CI_ub,[],tx_selected_table.gain(1))];
%         else
%             mRNA_CI_lb = [fluo_mRNA_conversion(tx_selected_table.CI_lb(1:new_split_point_idx),[],tx_selected_table.gain(1));...
%                 fluo_mRNA_conversion(tx_selected_table.CI_lb(new_split_point_idx+1:end),[],tx_selected_table.gain(new_split_point_idx + 1))];
%             mRNA_CI_ub = [fluo_mRNA_conversion(tx_selected_table.CI_ub(1:new_split_point_idx),[],tx_selected_table.gain(1));...
%                 fluo_mRNA_conversion(tx_selected_table.CI_ub(new_split_point_idx+1:end),[],tx_selected_table.gain(new_split_point_idx + 1))];
%         end
%         tx_selected_table = addvars(tx_selected_table,mRNA_CI_lb,'NewVariableNames','mRNA_concentration_CI_lb'); 
%         tx_selected_table = addvars(tx_selected_table,mRNA_CI_ub,'NewVariableNames','mRNA_concentration_CI_ub'); 
%     end
%     %%%%%% Protein %%%%%%
%     if ~isempty(PE_selected_table)
%         PE_split_point_idx_list = find(PE_selected_table.gain==PE_selected_table.gain(1)); 
%         PE_split_point_idx = PE_split_point_idx_list(end); 
%         if isequal(length(PE_split_point_idx_list),height(PE_selected_table))% if there's only one gain 
%             GFP_concentration = [fluo_GFP_conversion(PE_selected_table.fluorescence,[],PE_selected_table.gain(1))];
%         else
%                 % Applies to when there are 2 gains 
%             GFP_concentration = [fluo_GFP_conversion(PE_selected_table.fluorescence(1:PE_split_point_idx),[],PE_selected_table.gain(1));...
%                 fluo_GFP_conversion(PE_selected_table.fluorescence(PE_split_point_idx+1:end),[],PE_selected_table.gain(PE_split_point_idx + 1))];
%         end
%         PE_selected_table = addvars(PE_selected_table,GFP_concentration,'NewVariableNames','GFP_concentration'); 
%             % Replace negative values with nan 
%         PE_selected_table = remove_neg_val(PE_selected_table);
%             % Find new split after negative value removal
%         new_split_point_idx_list = find(PE_selected_table.gain==PE_selected_table.gain(1)); 
%         new_split_point_idx = new_split_point_idx_list(end); 
%         if isequal(length(new_split_point_idx_list),height(PE_selected_table))
%             GFP_CI_lb = [fluo_GFP_conversion(PE_selected_table.CI_lb,[],PE_selected_table.gain(1))];
%             GFP_CI_ub = [fluo_GFP_conversion(PE_selected_table.CI_ub,[],PE_selected_table.gain(1))];
%         else
%             GFP_CI_lb = [fluo_GFP_conversion(PE_selected_table.CI_lb(1:new_split_point_idx),[],PE_selected_table.gain(1));...
%                 fluo_GFP_conversion(PE_selected_table.CI_lb(new_split_point_idx+1:end),[],PE_selected_table.gain(new_split_point_idx + 1))];
%             GFP_CI_ub = [fluo_GFP_conversion(PE_selected_table.CI_ub(1:new_split_point_idx),[],PE_selected_table.gain(1));...
%                 fluo_GFP_conversion(PE_selected_table.CI_ub(new_split_point_idx+1:end),[],PE_selected_table.gain(new_split_point_idx + 1))];
%         end
%         PE_selected_table = addvars(PE_selected_table,GFP_CI_lb,'NewVariableNames','GFP_concentration_CI_lb'); 
%         PE_selected_table = addvars(PE_selected_table,GFP_CI_ub,'NewVariableNames','GFP_concentration_CI_ub'); 
%     end
%        % Consolidate columns/variables in PE_selected_table and
%        % tx_selected_table
%     tx_null_var = zeros(height(tx_selected_table),1);
%     PE_null_var = zeros(height(PE_selected_table),1); 
% 
%     tx_selected_table_vars = tx_selected_table.Properties.VariableNames; 
%     PE_selected_table_vars = PE_selected_table.Properties.VariableNames;
%     all_table_vars = [tx_selected_table_vars,PE_selected_table_vars];
%     all_table_vars = unique(all_table_vars);
% 
%     for var_idx = 1:length(all_table_vars)
%         var_name = all_table_vars{var_idx};
%         if ~any(strcmp(var_name,tx_selected_table_vars))
%             tx_selected_table  = addvars(tx_selected_table,tx_null_var,'NewVariableNames',var_name);
%         end
%         if ~any(strcmp(var_name,PE_selected_table_vars))
%             PE_selected_table  = addvars(PE_selected_table,PE_null_var,'NewVariableNames',var_name);
%         end
%     end
% 
%     % Reindex group number before combining 
%     num_group_tx_table = length(tx_data_description);
%     PE_selected_table.Group = PE_selected_table.Group + num_group_tx_table; 
%     selected_table = [tx_selected_table;PE_selected_table];
% 
%     % Add to problemObject
%     grouped_data = groupedData(selected_table); % define a groupedData object 
%     grouped_data.Properties.IndependentVariableName = 'Time'; 
%     grouped_data.Properties.GroupVariableName = 'Group'; 
% 
% end
%% Test adding parameters & rules to Simbiology model 
% test_species_name = get(test_Mobj.Species,'Name');
% RNA_species_idx = contains(test_species_name,'RNA utr') & ~contains(test_species_name,'RNase') & ~contains(test_species_name,'Ribo'); 
% RNA_species_name = test_species_name(RNA_species_idx); 
% test_reaction_name = get(test_Mobj.Reactions,'Reaction'); 
% tx_reaction_name_list = {};
% for test_rxn_name_idx = 1:length(test_reaction_name)
%     test_reaction_name_single = test_reaction_name{test_rxn_name_idx}; 
%     rhs_start_idx = find(test_reaction_name_single=='->');
%     truncated_test_reaction_name_single = test_reaction_name_single(rhs_start_idx:end); 
%     for RNA_species_idx = 1:length(RNA_species_name)
%         RNA_species_name_single = RNA_species_name{RNA_species_idx}; 
%         if contains(truncated_test_reaction_name_single,strcat('[',RNA_species_name_single,']'))
%             tx_reaction_name_list{end+1} = test_reaction_name_single;
%         end
%     end
% 
% end

%% Analyze initial results from adding toxin in fitting T7 strong no empty 
% result_file = load('param_est_run_save/20230914_param_est_run1610_.mat'); 
% SSE_summary = nan(48,1); 
% for iter = 1:48
%     fitResult = result_file.all_fitResults{iter,1}; 
%     % if ~isempty(fitResult)
%     SSE_summary(iter,1) = fitResult.SSE; 
%     % end
% end
% [~,min_idx] = min(SSE_summary); 
% opt_fitResult = result_file.all_fitResults{min_idx,1};
% params = [opt_fitResult.ParameterEstimates.Estimate];
% param_names = [opt_fitResult.ParameterEstimates.Name];
% 
% additional_track_species = {'toxin','AGTP'};
% simFunction = create_simFun_from_problemObject(result_file.problemObject,additional_track_species); 
% 
% tStart = 0; 
% tEnd = 21600;
% dosing_information = create_dosing_info_from_problemObject(result_file.problemObject);
% 
% 
% % Plot simulated data 
% modify_param_name_list = {'t7RNAP_0','AGTP_deg_time_0','tx_capacity_param'};
% modify_param_val_list = [10,1500,0.01]; 
% modified_params = modify_params(param_names,params,modify_param_name_list,modify_param_val_list); 
% 
% [Time,Data] = simFunction(modified_params',tEnd,dosing_information,tStart:tEnd);
% 
% % Plot simulated data 
% figure; 
% for data_idx = 1:length(Data)
%     timeVec = Time{data_idx,1}; 
%     dataVec = Data{data_idx,1}; 
%     plot(timeVec,dataVec(:,3),'LineWidth',1.5)
%     hold on 
% end


% simData = fitted(opt_fitResult); 
% 
% 
% figure; 
% for data_idx = 1:length(simData)
%     plot(simData(data_idx,1).Time,simData(data_idx,1).Data(:,24),'LineWidth',1.5)
%     hold on 
% end


% Check fitting for all 48 reusult structs
% observable_name = 'protein sfGFP*'; 
% for iter = 1:48
%     subplot(7,7,iter)
%     fitResult = result_file.all_fitResults{iter,1}; 
%     simData = fitted(fitResult); 
%     observable_idx = strcmp(simData(1,1).DataNames,observable_name);
%     for conc_idx = 1:length(conc_vec)
%         GFP_data = simData(conc_idx,1).Data(:,observable_idx); 
%         time_vec = simData(conc_idx,1).Time;
%         plot(time_vec,GFP_data,'LineWidth',1.5)
%         hold on
%     end
% end

% Toxin time-course
% figure;
% for simData_idx = 1:length(Data)
%     subplot(3,3,simData_idx)
%     toxin_time_course = simData(simData_idx,1).Data(:,99);
%     [toxin_time_course_unique,ia,ic] = unique(toxin_time_course); 
%     time_course = simData(simData_idx,1).Time(ia); 
%     even_timeVec = 0:10:time_course(1); 
%     interp_toxin_time_course = interp1(time_course,toxin_time_course_unique,even_timeVec','linear','extrap'); 
%     plot(simData(simData_idx,1).Time,simData(simData_idx,1).Data(:,99),'LineWidth',1.5)
%     hold on 
% end

%% Plot PE baseline fitting results 
% result_analysis_file = load('param_est_run_save/20230813_param_est_run0807_PE_sigma70_weak_no_empty.mat'); 
% [~,opt_run_idx] = min(result_analysis_file.metric_summary); 
% fit_result_file = load(sprintf('param_est_run_save/20230813_param_est_run0807_PE_sigma70_weak_no_empty_fitResult%d.mat',opt_run_idx)); 
% simData = fitted(fit_result_file.fitResults); 
% 
%% sigma70 weak resource check 
% resource_name_list = {'t7RNAP','RNAP','RNase','AGTP','RNA'};
% data_name_list = {'RNAP','RNase','Ribo','AGTP','RNA utrGFP--sfGFP','RNAP:DNA pJ23105--utrGFP--sfGFP','RNAP:DNA pJ23119--utrGFP--sfGFP',...
%     'RNAP:DNA pkanR--utrkanR--kanR','RNA utrkanR--kanR','protein kanR','protein sfGFP*'};

% for resource_idx = 1:length(resource_name_list)
%     resource_name = resource_name_list{resource_idx}; 
%     data_idx_list = get_data_idxlist(simData(15).DataNames,resource_name); 
% data_idx_list = nan(length(data_name_list),1); 
% for i = 1:length(data_name_list)
%     data_idx_list(i) = find(strcmp(simData(1).DataNames,data_name_list{i}));
% 
% end

% % Increase RNAP concentration to alleviate competition
% [~,min_idx] = min(result_analysis_file.metric_summary);
% opt_fitResults = result_analysis_file.all_fitResults{min_idx,1}; 
% 
% param_names = opt_fitResults.ParameterEstimates.Name; 
% param_vals = opt_fitResults.ParameterEstimates.Estimate; 
% 
% RNAP_idx = find(strcmp(param_names,'RNAP_0'));
% param_vals(RNAP_idx) = 30; 
% 
% Mobj = result_analysis_file.problemObject.Model; 
% fitted_Mobj = assign_model_parameters(Mobj,param_names,param_vals); 
% 
% model_plasmid_name_list = {'DNA pJ23105--utrGFP--sfGFP','DNA pkanR--utrkanR--kanR'};
% simData = cell(7,1); 
% for conc_idx = 1:length(conc_vec)
%     conc = conc_vec(conc_idx); 
%     temp_model = copyobj(fitted_Mobj);
%     for plasmid_idx = 1:length(model_plasmid_name_list)
%         model_plasmid_name = model_plasmid_name_list{plasmid_idx}; 
%         species_obj = sbioselect(temp_model,'Type','Species','Name',model_plasmid_name);
%         set(species_obj,'Value',conc);
%     end
%     simData{conc_idx,1} = txtl_runsim(temp_model,14*60*60); 
% end


% Check whether parameters in reaction -> KineticLaw object are all deleted
% They are .
% test_model = fit_result_file.problemObject.Model; 
% Mobj_rxns = get(test_model,'Reactions');
% for rxn_idx = 1:length(Mobj_rxns)
%     rxn = Mobj_rxns(rxn_idx);
%     rxn_params = rxn.KineticLaw.Parameters;
%     for rxn_param_idx = 1:length(rxn_params)
%         rxn_param = rxn_params(rxn_param_idx); 
%         model_param = sbioselect(test_model.Parameters,'Name',rxn_param.Name); 
%         if ~isemtpy(model_param)
%             fprintf(rxn_param.Name)
%             fprintf(rxn_param.Value == model_param.Value)
%         end
%     end
% end

% figure; 
% for plot_data_name_idx = 1:length(data_idx_list)
%     subplot_num = ceil(sqrt(length(data_idx_list)));
%     subplot(subplot_num,subplot_num,plot_data_name_idx)
%     for conc_idx = 1:7
%         single_simData = simData(conc_idx); 
%         if max(single_simData.Data(:,data_idx_list(plot_data_name_idx))) > 1e-08
%             plot(single_simData.Time,single_simData.Data(:,data_idx_list(plot_data_name_idx)),'LineWidth',1.5)%,'Color','g')
%             hold on 
%         end
%         title(single_simData.DataNames{data_idx_list(plot_data_name_idx)})
%     end
% end
% legend({'0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM'})
% sgtitle('sigma70 weak')


% Plot simulated data 
% figure;
% for prom_idx = 1:1
%     subplot(2,2,prom_idx)
%     for conc_idx = 1:7
%         % single_simData = simData((prom_idx - 1) * 7 + conc_idx);
%         single_simData = simData{conc_idx};
%         data_idx = find(strcmp(single_simData.DataNames,'protein sfGFP*')); 
%         plot(single_simData.Time,single_simData.Data(:,data_idx),'LineWidth',1.5)
%         hold on 
%         % title(promotors{prom_idx})
%         title('sigma70 weak')
%     end
% end
% legend({'0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM'})
% sgtitle('Simulated')

% Plot experimental data 
% experimental_data = fit_result_file.problemObject.Data; 
% 
% figure; 
% for prom_idx = 1:4
%     subplot(2,2,prom_idx)
%     for conc_idx = 1:7
%         selected_table = experimental_data(experimental_data.Group == 84 + (prom_idx - 1) * 28 + conc_idx,:); 
%         plot(selected_table.Time,selected_table.GFP_concentration,'LineWidth',1.5)
%         hold on 
%         title(promotors{prom_idx})
%     end
% end
% legend({'0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM'})
% sgtitle('Experimental')
%% Investigate lack of huge positive crosstalk in sigma70 weak 
% crosstalk_trained_file = load('param_est_run_save/20230808_param_est_run0910_.mat'); 
% all_fitResults = crosstalk_trained_file.all_fitResults; 
% simData = fitted(all_fitResults{17,1});
% start_dataset_idx = 85; 
% plot_all_species_time_course(simData,start_dataset_idx)
% 
% crosstalk_ratio = calculate_crosstalk_ratio(simData,7,4,'PE'); 
% figure; 
% for promotor_idx = 1:4
%     subplot(2,2,promotor_idx)
%     plot(conc_vec,crosstalk_ratio{1,promotor_idx},'Color',[0.5,0.5,0.5]);
%     hold on 
%     plot(conc_vec,crosstalk_ratio{2,promotor_idx},'Color','r');
%     plot(conc_vec,crosstalk_ratio{3,promotor_idx},'Color','b');
% 
% end
% sgtitle('crosstalk fitting fitting #17')


%% Check error model & short time-course fitting results 
% modes = {'TX','PE','combined'};
% % 
% % for mode_idx = 3:length(modes)
% %     mode = modes{mode_idx}; 
% %     all_fitResults = cell(48,1); 
% %     for rep = 1:48
% %         % try
% %         % fitResult_file = load(sprintf('param_est_run_save/20230814_param_est_run1720_%s_error_model_fitResult%d.mat',mode,rep));
% %         % fit_result_path = sprintf('param_est_run_save/20230814_param_est_run1720_%s',mode); 
% % 
% %         fitResult_file = load(sprintf('param_est_run_save/20230814_param_est_run1221_%s_short_time_course_fitResult%d.mat',mode,rep));
% %         fit_result_path = sprintf('param_est_run_save/20230814_param_est_run1221_short_time_course_%s',mode); 
% %         % catch
% %         %     fitResult_file = load(sprintf('param_est_run_save/20230813_param_est_run0751_%s_fitResult%d.mat',label,rep));
% %         %     fit_result_path = sprintf('param_est_run_save/20230813_param_est_run0751_%s',label); 
% %         % end
% %         all_fitResults{rep,1} = fitResult_file.fitResults;
% %         if isequal(rep,1)
% %             problemObject = fitResult_file.problemObject; 
% %         end
% %     end
% % 
% %     metric_name = 'SSE'; 
% %     plot_fitting = true;
% %     plot_param_dist = true;
% %     wrapper_analyze_fit_result_v2(all_fitResults,problemObject,fit_result_path,metric_name,plot_fitting,plot_param_dist);
% % end
% for mode_idx = 2:2
%     mode = modes{mode_idx}; 
%     fit_result_analysis_path = sprintf('param_est_run_save/20230814_param_est_run1221_short_time_course_%s.mat',mode); 
%     fit_result_analysis_file = load(fit_result_analysis_path); 
% 
%     [~,opt_idx] = min(fit_result_analysis_file.metric_summary); 
%     fit_result_object = load(sprintf('param_est_run_save/20230814_param_est_run1221_%s_short_time_course_fitResult%d.mat',mode,opt_idx));
%     simulated_data = fitted(fit_result_object.fitResults); 
%     num_conc = length(conc_vec); 
%     num_promotor = 4; 
% 
%     crosstalk_ratio = calculate_crosstalk_ratio(simulated_data,num_conc,num_promotor,mode); 
% 
%     figure; 
%     for promotor_idx = 1:num_promot
%         subplot(2,2,promotor_idx)
%         plot(conc_vec,crosstalk_ratio{1,promotor_idx},'Color',[0.5,0.5,0.5]);
%         hold on 
%         plot(conc_vec,crosstalk_ratio{2,promotor_idx},'Color','r');
%         plot(conc_vec,crosstalk_ratio{3,promotor_idx},'Color','b');
% 
%     end
%     sgtitle(sprintf('%s crosstalk ratio weight %d',mode,weight))
%     % saveas(gcf,sprintf('plots/20230815_test_weight_%s_crosstalk_ratio_weight %d',mode,weight),'png'); 
% 
% end
%% Check weights on t = 3hrs (protein) or t_maxval (mRNA)
% all_weights = [1e+03,1e+04,1e+05,1e+06,1e+07,1e+08];
% all_modes = {'TX','PE','combined'};
% num_conc = 7;
% num_promotor = 4; 
% conc_vec = [0.5,1,2.5,5,10,15,30]; 

% Check fitting 
% for mode_idx = 1:length(all_modes)
%     mode = all_modes{mode_idx}; 
%     for weight_idx = 1:length(all_weights)
%         weight = all_weights(weight_idx);
%         all_fitResults = cell(48,1); 
%         for rep = 1:48
%             try
%                 fitResult_file = load(sprintf('param_est_run_save/20230814_param_est_run1222_%s_%d_fitResult%d.mat',mode,weight,rep)); 
%                 fit_result_path = sprintf('param_est_run_save/20230814_param_est_run1222_%s_%d.mat',mode,weight);
%             catch
%                 fitResult_file = load(sprintf('param_est_run_save/20230814_param_est_run1223_%s_%d_fitResult%d.mat',mode,weight,rep)); 
%                 fit_result_path = sprintf('param_est_run_save/20230814_param_est_run1223_%s_%d.mat',mode,weight);
%             end
%             all_fitResults{rep,1} = fitResult_file.fitResults;
%             if isequal(rep,1)
%                 problemObject = fitResult_file.problemObject; 
%             end
%         end
%         metric_name = 'SSE'; 
%         plot_fitting = true;
%         plot_param_dist = true;
%         wrapper_analyze_fit_result_v2(all_fitResults,problemObject,fit_result_path,metric_name,plot_fitting,plot_param_dist); 
% 
% 
%     end
% end

% % Check crossfit ratio 
% for mode_idx = 1:3
%     mode = all_modes{mode_idx}; 
%     for weight_idx = 1:length(all_weights)
%         weight = all_weights(weight_idx);
%         try
%             fit_analysis_fileName = sprintf('param_est_run_save/20230814_param_est_run1222_%s_%d.mat',mode,weight);
%             fit_analysis_file = load(fit_analysis_fileName);
%         catch
%             fit_analysis_fileName = sprintf('param_est_run_save/20230814_param_est_run1223_%s_%d.mat',mode,weight); 
%             fit_analysis_file = load(fit_analysis_fileName);
%         end
%         metric_summary = fit_analysis_file.metric_summary; 
% 
%         % Default to check best fit 
%         [~,min_idx] = min(metric_summary); 
%         opt_fit_result_file = load(sprintf('%s_fitResult%d.mat',fit_analysis_fileName(1:end-4),min_idx)); 
%         simulated_data = fitted(opt_fit_result_file.fitResults); 
%         if ~strcmp(mode,'combined')
%             if strcmp(mode,'TX')
%                 num_promotor = 3;
%             elseif strcmp(mode,'PE')
%                 num_promotor = 4;
%             end
%             crosstalk_ratio = calculate_crosstalk_ratio(simulated_data,num_conc,num_promotor,mode); 
%             figure; 
%             for promotor_idx = 1:num_promotor
%                 subplot(2,2,promotor_idx)
%                 plot(conc_vec,crosstalk_ratio{1,promotor_idx},'Color',[0.5,0.5,0.5]);
%                 hold on 
%                 plot(conc_vec,crosstalk_ratio{2,promotor_idx},'Color','r');
%                 plot(conc_vec,crosstalk_ratio{3,promotor_idx},'Color','b');
% 
%             end
%             sgtitle(sprintf('%s crosstalk ratio weight %d',mode,weight))
%             saveas(gcf,sprintf('plots/20230815_test_weight_%s_crosstalk_ratio_weight %d',mode,weight),'png'); 
%         else
%             TX_crosstalk_ratio = calculate_crosstalk_ratio(simulated_data,num_conc,3,'TX'); 
%             PE_crosstalk_ratio = calculate_crosstalk_ratio(simulated_data,num_conc,4,'PE'); 
%             figure; 
%             for promotor_idx = 1:3
%                 subplot(2,2,promotor_idx)
%                 plot(conc_vec,TX_crosstalk_ratio{1,promotor_idx},'Color',[0.5,0.5,0.5]);
%                 hold on 
%                 plot(conc_vec,TX_crosstalk_ratio{2,promotor_idx},'Color','r');
%                 plot(conc_vec,TX_crosstalk_ratio{3,promotor_idx},'Color','b');
% 
%             end
%             sgtitle(sprintf('Combined mode TX crosstalk ratio weight %d',weight))
%             saveas(gcf,sprintf('plots/20230815_test_weight_combined_TX_crosstalk_ratio_weight %d',weight),'png'); 
%             figure; 
%             for promotor_idx = 1:4
%                 subplot(2,2,promotor_idx)
%                 plot(conc_vec,crosstalk_ratio{1,promotor_idx},'Color',[0.5,0.5,0.5]);
%                 hold on 
%                 plot(conc_vec,crosstalk_ratio{2,promotor_idx},'Color','r');
%                 plot(conc_vec,crosstalk_ratio{3,promotor_idx},'Color','b');
% 
%             end
%             sgtitle(sprintf('Combined mode PE crosstalk ratio weight %d',weight))
%             saveas(gcf,sprintf('plots/20230815_test_weight_combined_PE_crosstalk_ratio_weight %d',weight),'png'); 
%         end
%     end
% end
%% Check subset fitting results 

% subset_info_file = load('subset_fitting_descriptions.mat'); 
% labels = subset_info_file.labels; 
% 
% 
% % timepoint_choice = {'0751','0807'};
% for subset_idx = 13:length(labels)
%     label = labels{subset_idx}; 
%     all_fitResults = cell(48,1); 
%     for rep = 1:48
%         try
%             fitResult_file = load(sprintf('param_est_run_save/20230813_param_est_run0807_%s_fitResult%d.mat',label,rep));
%             fit_result_path = sprintf('param_est_run_save/20230813_param_est_run0807_%s',label); 
%         catch
%             fitResult_file = load(sprintf('param_est_run_save/20230813_param_est_run0751_%s_fitResult%d.mat',label,rep));
%             fit_result_path = sprintf('param_est_run_save/20230813_param_est_run0751_%s',label); 
%         end
%         all_fitResults{rep,1} = fitResult_file.fitResults;
%         if isequal(rep,1)
%             problemObject = fitResult_file.problemObject; 
%         end
%     end
% 
%     metric_name = 'SSE'; 
%     plot_fitting = true;
%     plot_param_dist = true;
%     wrapper_analyze_fit_result_v2(all_fitResults,problemObject,fit_result_path,metric_name,plot_fitting,plot_param_dist); 
% 
% end



%% Generate subsets of group description & save somewhere 

% base_group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%         'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%         'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70',...
%         'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%         'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%         'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%         'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% 
% all_group_descriptions = {}; 
% labels = {}; 
% 
% % Just singles for the 1st round 
% for i = 1:length(base_group_description)
%     all_group_descriptions = [all_group_descriptions {base_group_description(i)}];
%     labels = [labels base_group_description(i)];
% 
% end
% % Baselines for the 2nd round 
% TX_baseline_description = base_group_description(contains(base_group_description,'no_empty') & ~contains(base_group_description,'PE')); 
% PE_baseline_description = base_group_description(contains(base_group_description,'no_empty') & contains(base_group_description,'PE')); 
% combined_baseline_description = [TX_baseline_description,PE_baseline_description];
% all_group_descriptions = [all_group_descriptions {TX_baseline_description} {PE_baseline_description} {combined_baseline_description}];
% labels = [labels {'TX_baseline','PE_baseline','combined_baseline'}]; 
% 
% % Promotor specific for the 3rd round 
% promotor_list = {'T7_strong','T7_weak','sigma70_strong','sigma70_weak'};
% for promotor_idx = 1:length(promotor_list)
%     promotor_name = promotor_list{promotor_idx}; 
%     TX_description = base_group_description(contains(base_group_description,promotor_name) & ~contains(base_group_description,'PE'));
%     PE_description = base_group_description(contains(base_group_description,promotor_name) & contains(base_group_description,'PE'));
%     combined_descrption = [TX_description,PE_description];
%     if isequal(promotor_idx,length(promotor_list))
%         all_group_descriptions = [all_group_descriptions {PE_description}];
%         labels = [labels {sprintf('%s_PE',promotor_name)}];
%     else
%         all_group_descriptions = [all_group_descriptions {TX_description} {PE_description} {combined_descrption}];
%         labels = [labels {sprintf('%s_TX',promotor_name),sprintf('%s_PE',promotor_name),sprintf('%s_combined',promotor_name)}];
%     end
% 
% end
% 
% save('subset_fitting_descriptions','all_group_descriptions','labels');


%% Resource utilization and 
% load('param_est_run_save/20230809_param_est_run1440_.mat'); 

% min_SSE_run42 = all_fitResults{42,1};
% large_SSE_run46 = all_fitResults{46,1}; 
% fitted_data_42 = fitted(min_SSE_run42);
% fitted_data_46 = fitted(large_SSE_run46); 
% 
% % simData #85 - 0.5nM sigma70 weak no empty 
% % simData #92 - 0.5nM sigma70 weak empty 
% % simData #99 - 0.5nM sigma70 weak empty T7
% % simData #106 - 0.5nM sigma70 weak empty sigma70
%     % Plot these on the same plot in different colors 
%     % Input: start_dataset_idx
% start_dataset_idx = 89; 
% plot_all_species_time_course(fitted_data_46,start_dataset_idx); 
%% Check OptimResult
% param_est_result_file_w_broc_param = load('param_est_run_save/20230802_param_est_run1352_.mat'); 
% param_est_result_file_wo_broc_param = load('param_est_run_save/20230808_param_est_run0910_.mat'); 
% 
% 
% w_metric_summary = param_est_result_file_w_broc_param.metric_summary;
% wo_metric_summary = param_est_result_file_wo_broc_param.metric_summary;
% w_param_summary = param_est_result_file_w_broc_param.estimated_param_summary;
% wo_param_summary = param_est_result_file_wo_broc_param.estimated_param_summary;
% 
% for rep = 1:48
%     figure; 
%     SSE_w = w_metric_summary(rep); 
%     SSE_wo = wo_metric_summary(rep); 
%     for param_idx = 1:height(w_param_summary)
%         subplot(6,6,param_idx)
% 
% 
%         w_param_name = w_param_summary.Name(param_idx);
%         wo_param_idx = strcmp(w_param_name,wo_param_summary.Name);
% 
%         eval(sprintf('w_paramVal = w_param_summary.Estimate_%d(param_idx);',rep));
%         eval(sprintf('wo_paramVal = wo_param_summary.Estimate_%d(wo_param_idx);',rep));
% 
%         bar([w_paramVal,wo_paramVal]); 
%         title(sprintf('%s',strrep(w_param_name{1,1},'_',' '))); 
% 
% 
%     end
%     sgtitle(sprintf('SSE_w = %d, SSE_wo = %d',SSE_w,SSE_wo))
% end

%% Load result file 
% baseline_trained_file = load('param_est_run_save/20230802_param_est_run1459_.mat');
% crosstalk_trained_file = load('param_est_run_save/20230802_param_est_run1352_.mat'); 
% tx_baseline_trained_file = load('param_est_run_save/20230509_param_est_run1334_.mat');
% tx_crosstalk_trained_file = load('param_est_run_save/20230509_param_est_run1509_.mat'); 
% % combined_fitResult = cell(48,1); 
% % for i = 1:48
% %     combined_data_trained_fitResult_file = load(sprintf('param_est_run_save/20230804_multiple_fitting_fitResult%d.mat',i)); 
% %     combined_fitResult{i,1} = combined_data_trained_fitResult_file.fitResult; 
% % end
% combined_result_analysis = load('param_est_run_save/20230804_param_est_run1413_.mat'); 

%% Compare best-fitting parameters from baseline and crosstalk (and the tx-level ones) up
% [PE_baseline_param_names,PE_baseline_opt_param] = get_best_fitted_param(baseline_trained_file); 
% [PE_crosstalk_param_names,PE_crosstalk_opt_param] = get_best_fitted_param(crosstalk_trained_file); 
% [tx_baseline_param_names,tx_baseline_opt_param] = get_best_fitted_param(tx_baseline_trained_file); 
% [tx_crosstalk_param_names,tx_crosstalk_opt_param] = get_best_fitted_param(tx_crosstalk_trained_file); 
% [combined_crosstalk_param_names,combined_crosstalk_opt_param] = get_best_fitted_param(combined_result_analysis); 
% 
%     % Use PE crosstalk as it should have the most fitted parameters 
% X_labels = {'Broc baseline','Broc crosstalk','PE baseline','PE crosstalk','Combined crosstalk'};
% for param_idx = 1:length(combined_crosstalk_param_names)
%     % 
%     subplot(ceil(sqrt(length(combined_crosstalk_param_names))),ceil(sqrt(length(combined_crosstalk_param_names))),param_idx)
% 
%     % Get parameter value for each fitted result 
%     param_name = combined_crosstalk_param_names{param_idx}; 
%     combined_crosstalk_paramVal = combined_crosstalk_opt_param(param_idx); 
% 
%     PE_crosstalk_paramVal = PE_crosstalk_opt_param(strcmp(param_name,PE_crosstalk_param_names)); 
%     PE_baseline_paramVal = PE_baseline_opt_param(strcmp(param_name,PE_baseline_param_names)); 
%     tx_baseline_paramVal = tx_baseline_opt_param(strcmp(param_name,tx_baseline_param_names));
%     tx_crosstalk_paramVal = tx_crosstalk_opt_param(strcmp(param_name,tx_crosstalk_param_names)); 
% 
%     % Select x-axis labels 
%     param_array = {tx_baseline_paramVal,tx_crosstalk_paramVal,PE_baseline_paramVal,PE_crosstalk_paramVal,combined_crosstalk_paramVal}; 
%     temp_X_labels = {}; 
%     for empty_idx = 1:length(param_array)
%         if ~isempty(param_array{1,empty_idx})
%             temp_X_labels = [temp_X_labels X_labels{empty_idx}];
%         end
%     end
% 
%     X = categorical(temp_X_labels);
%     X = reordercats(X,temp_X_labels); 
% 
%     bar(X,[tx_baseline_paramVal,tx_crosstalk_paramVal,PE_baseline_paramVal,PE_crosstalk_paramVal,combined_crosstalk_paramVal]); 
% 
%     title(strrep(param_name,'_',' '));    
% end
% sgtitle('Optimal parameter value estimated with different training data'); 


%% Check non-optimal fitted models and their single-value crosstalk ratio 
    % Given two result files, compare the two parameter tables 
    % Plot out each parameter value with their corresponding 3 crosstalk
    % ratios 

    % Figure 1: Plot the line for crosstalk ratio for each fitResult 
% crosstalk_trained_file = load('param_est_run_save/20230802_param_est_run1352_.mat'); 
% crosstalk_trained_file = load('param_est_run_save/20230809_param_est_run1440_.mat'); 
%     % Get GFP concentration at t = 10800s 
% num_conc = 7;
% conc_vec = [0.5,1,2.5,5,10,15,30]; 
% num_promotor = 4;
% for i = 1:length(crosstalk_trained_file.all_fitResults)
%     fit_result = crosstalk_trained_file.all_fitResults{i}; 
%     fitted_param_summary = crosstalk_trained_file.estimated_param_summary; 
%     num_fitted_parameters = height(fitted_param_summary); 
% 
%     simulated_data = fitted(fit_result); 
%     crosstalk_ratio = calculate_crosstalk_ratio(simulated_data,num_conc,num_promotor,'PE'); 
% 
%     figure; 
%     for promotor_idx = 1:num_promotor
%         subplot(2,2,promotor_idx)
%         plot(conc_vec,crosstalk_ratio{1,promotor_idx},'Color',[0.5,0.5,0.5]);
%         hold on 
%         plot(conc_vec,crosstalk_ratio{2,promotor_idx},'Color','r');
%         plot(conc_vec,crosstalk_ratio{3,promotor_idx},'Color','b');
% 
%     end
%     sgtitle(sprintf('Crosstalk ratio fitting #%d',i))
% 
% end

function [param_names,opt_val] = get_best_fitted_param(result_file)
    [~,opt_idx] = min(result_file.metric_summary); 
    param_summary = result_file.estimated_param_summary; 
    param_names = param_summary.Name; 
    eval(sprintf('opt_val = param_summary.Estimate_%d;',opt_idx)); 
    
end


%% Troubleshoot zero TL_translation_kanR mystery
    % Issue: This parameter should be assigned in model.rules. The relevant
    % rule shows up after running setProblemObject_updated and fit, but not
    % after runnning setProblemObject_PE and fit. Let's do these things in
    % parallel and see what's wrong 
% group_description = {'sigma70_strong_no_empty'}; % Give a small group of data
%     % I don't think this matters as I swapped this out and didn't seem to
%     % impact anything 
% param_info_path = '/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/benchmark_model/parameters_test2.xlsx';
% 
%     % Get the problemObject 
% problemObject_updated = setProblemObject_updated(group_description,param_info_path); 
% problemObject_PE = setProblemObject_PE(group_description,param_info_path); 
% 
%     % fit the problemObject
% for iter = 1:1
%     % Modify initial conditions, sample within bounds
%     rng(iter)
%     temp_problemObject_updated = problemObject_updated; 
%     temp_problemObject_PE = problemObject_PE; 
%     [fitResults_updated,simdataI_updated] = fit(temp_problemObject_updated);
%     [fitResults_PE,simdataI_PE] = fit(temp_problemObject_PE);
% end
% simData_updated = fitted(fitResults_updated);
% simData_PE = fitted(fitResults_PE); 
% 
% 
% test_fit_results2 = load('param_est_run_save/20230509_param_est_run1334_.mat'); 
%% Add a variable to mRNA data table with sigma70 weak plasmid concentration
    % TO make mRNA data table have the same variables as PE data table 
% 
% tx_data_file = load('data_structures/simbio_data_table_updated_FP.mat'); 
% variable_name_stem_list = tx_data_file.variable_name_stem_list;
% data_table = tx_data_file.all_data_table; 
% sigma70_weak_reporter_plasmid = zeros(height(data_table),1); 
% all_data_table = addvars(data_table,sigma70_weak_reporter_plasmid,'NewVariableNames',...
%     'sigma70_weak_reporter_plasmid','After','sigma70_strong_reporter_plasmid'); 
% save('data_structures/simbio_data_table_updated_FP.mat','all_data_table','variable_name_stem_list');
%% Plot Figure 2 in plasmid crosstalk paper from processed data 
% promotor_list = {'T7_strong','T7_weak','sigma70_strong'};
% promotor_list_old = {'T7_strong','T7_weak','sigma70'};
% empty_list = {'no_empty','empty','empty_T7','empty_sigma70'}; 
% color_list = {'g',[0.5,0.5,0.5],'r','b'};
% appendix_list = {'','_updated_FP','_updated_AP'}; 
% style_list = {':','-','-'}; 
% conc_list = [0.5,1,2.5,5,10,15,30];
% 
% % appendix = appendix_list{3}; 
% 
% for appendix_idx = 2:length(appendix_list)
%     appendix = appendix_list{appendix_idx};
%     figure; 
%     for conc_idx = 1:length(conc_list)
%         conc = conc_list(conc_idx); 
%         for prom_idx = 1:length(promotor_list)
%             promotor_name = promotor_list{prom_idx}; 
%             subplot(length(conc_list),length(promotor_list),length(promotor_list) * (conc_idx - 1) + prom_idx)
%             hold on 
%             if isequal(conc_idx,1)
%                 title(strrep(promotor_name,'_',' '))
%                 hold on 
%             end
%             if isequal(prom_idx,1)
%                 ylabel(sprintf('%.2f nM',conc))
%                 hold on 
%             end
%             if isequal(conc_idx,length(conc_list))
%                 xlabel('time(s)')
%                 hold on 
%             end
%             for empty_idx = 1:length(empty_list)
%                 empty_name = empty_list{empty_idx}; 
%                 data_fileName = sprintf('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/processed_data/%s_%s%s.xlsx',promotor_name,empty_name,appendix);
%                 timeVec_cell = readcell(data_fileName,'Sheet','timeVec'); 
%                 fluorescence_data_cell = readcell(data_fileName,'Sheet','fluorescence_data'); 
%                 timeVec = cell2mat(timeVec_cell(2:end,2)); 
%                 fluorescence_data = cell2mat(fluorescence_data_cell(2:end,conc_idx + 1));
% 
%                 plot(timeVec,fluorescence_data,'Color',color_list{empty_idx},'LineWidth',1.5,'LineStyle',style_list{appendix_idx}); 
%                 hold on 
%             end
%         end
%     end
% end

% for dum_idx = 1:4
%     plot(nan,nan,'Color','g','LineWidth',1.5)
%     hold on 
%     plot(nan,nan,'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(nan,nan,'Color','r','LineWidth',1.5)
%     plot(nan,nan,'Color','b','LineWidth',1.5)
%     plot(nan,nan,'Color','g','LineWidth',1.5,'LineStyle','--')
%     plot(nan,nan,'Color',[0.5,0.5,0.5],'LineWidth',1.5,'LineStyle','--')
%     plot(nan,nan,'Color','r','LineWidth',1.5,'LineStyle','--')
%     plot(nan,nan,'Color','b','LineWidth',1.5,'LineStyle','--')
% end
% legend('No Empty FP','Empty FP','Empty T7 FP','Empty sigma70 FP','No Empty AP','Empty AP','Empty T7 AP','Empty sigma70 AP')

%%
% out = test_vararg(); 
% [out1,out2] = test_vararg();
% 
% function varargout = test_vararg(varargin)
%     varargout{1} = [];
%     varargout{2} = 'dfdfd'; 
% 
% end