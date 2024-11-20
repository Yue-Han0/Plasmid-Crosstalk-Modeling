clear
clc

currentpath = pwd; 
addpath(genpath(currentpath))
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))
% addpath('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling-/processed_data')
promotor_name_list = {'T7_strong','T7_weak','sigma70_strong','sigma70_weak'};
conc_vec = [0.5,1,2.5,5,10,15,30]; 
mode = 'PE';
num_iter = 48; % Unless otherwise specified, 48 initial conditions were run 
num_conc = 7;
num_promotor = 4; 
all_obj_string_file = load('qualitative_trend_obj_string.mat'); 
all_obj_string = all_obj_string_file.all_obj_string; 
tStart = 0; 
tEnd = 21600; % Modify this later to draw info from data  
reporter_species_name = 'protein sfGFP*'; 


%% Visualization of # samples that meet criteria

% init_param_file = load('test_save_files/20240819_lhs_sampled_params.mat'); 
% all_satisfied_sampled_params = init_param_file.satisfied_sampled_params; 
% sampling_result_file2 = load('test_save_files/202409_high_res_penalty_sampling_result_summary.mat'); 
% sample_sampling_resultFile2 = load('param_est_run_save/20240819_param_sampling_constrained_run10.mat');
% problemObject_sampling2 = sample_sampling_resultFile2.problemObject; 
% dosing_information_sampling2 = init_param_file.dosing_information; 
% 
%     % Create simFunction using problemObject
% all_model_species2 = {problemObject_sampling2.Model.Species.Name}; 
% reporter_species_idx2 = find(strcmp(all_model_species2,reporter_species_name)); 
% additional_track_species2 = [all_model_species2(1:reporter_species_idx2 - 1),all_model_species2(reporter_species_idx2 + 1:end)];
% simFunction_sampling2 = create_simFun_from_problemObject(problemObject_sampling2,additional_track_species2); 
% 
%     % Load penalty cutoff 
% high_res_penalty_goals = load('test_save_files/high_resolution_penalty_goals.mat'); 
% penalty_cutoff = high_res_penalty_goals.penalty_goals; 
% 
%     % Use plot_simulated to compare
% num_case = 5; 
% 
% %     % Calculate number of criteria each row meets
% % satisfy_flag_matrix = sampling_result_file2.all_high_res_penalty_terms < penalty_cutoff; 
% % num_criteria_satisfied = sum(satisfy_flag_matrix,2); 
% % milestones = [5,10,15,20,25,30,32];
% % num_satisfied_milestones = nan(1,length(milestones)); 
% % for milestone_idx = 1:length(milestones)
% %     milestone = milestones(milestone_idx); 
% %     num_satisfied_milestones(milestone_idx) = sum(num_criteria_satisfied >= milestone); 
% % end
% % labels = {'>= 5','>= 10','>= 15','>= 20','>= 25','>= 30','>= 32'};
% % funnelchart(num_satisfied_milestones,'Labels',labels)
% 
% % Sankey diagram for number of satisfied qualitative trends 
% 
%     % First define 9 criteria used
%         % (1) Baseline dev 
%         % (2) Crosstalk ratio dev
%         % (3) Baseline penalty (4 promotor cases)
%         % (4) T7 strong crosstalk (3 combo cases) 
%         % (5) T7 weak crosstalk (3 combo cases) 
%         % (6) sigma70 strong crosstalk (3 combo cases) 
%         % (7) sigma70 weak crosstalk (3 combo cases)
%         % (8) Larger positive crosstalk in sigma70 weak
%         % (9) Residual mRNA concentration 
% baseline_dev = sampling_result_file2.all_high_res_penalty_terms(:,1); 
% crosstalk_ratio_dev = sampling_result_file2.all_high_res_penalty_terms(:,2); 
% baseline_penalty = sampling_result_file2.all_high_res_penalty_terms(:,3:6); 
% 
% T7_strong_positive_crosstalk_penalty_idx_list = 7:4:15;
% T7_strong_negative_crosstalk_penalty_idx_list = 19:4:27;
% T7_weak_positive_crosstalk_penalty_idx_list = 8:4:16;
% T7_weak_negative_crosstalk_penalty_idx_list = 20:4:28;
% sigma70_strong_positive_crosstalk_penalty_idx_list = 9:4:17;
% sigma70_strong_negative_crosstalk_penalty_idx_list = 21:4:29;
% sigma70_weak_positive_crosstalk_penalty_idx_list = 10:4:18;
% sigma70_weak_negative_crosstalk_penalty_idx_list = 22:4:30;
% T7_strong_crosstalk_penalty = sampling_result_file2.all_high_res_penalty_terms(:,[T7_strong_positive_crosstalk_penalty_idx_list,T7_strong_negative_crosstalk_penalty_idx_list]); 
% T7_weak_crosstalk_penalty = sampling_result_file2.all_high_res_penalty_terms(:,[T7_weak_positive_crosstalk_penalty_idx_list,T7_weak_negative_crosstalk_penalty_idx_list]);     
% sigma70_strong_crosstalk_penalty = sampling_result_file2.all_high_res_penalty_terms(:,[sigma70_strong_positive_crosstalk_penalty_idx_list,sigma70_strong_negative_crosstalk_penalty_idx_list]); 
% sigma70_weak_crosstalk_penalty = sampling_result_file2.all_high_res_penalty_terms(:,[sigma70_weak_positive_crosstalk_penalty_idx_list,sigma70_weak_negative_crosstalk_penalty_idx_list]); 
% 
% large_positive_crosstalk_penalty = sampling_result_file2.all_high_res_penalty_terms(:,end-1);
% residual_mRNA_penalty = sampling_result_file2.all_high_res_penalty_terms(:,end); 
% 
% %     % Reorganize penalty cutoff
% T7_strong_crosstalk_penalty_cutoff = penalty_cutoff([T7_strong_positive_crosstalk_penalty_idx_list,T7_strong_negative_crosstalk_penalty_idx_list]); 
% T7_weak_crosstalk_penalty_cutoff = penalty_cutoff([T7_weak_positive_crosstalk_penalty_idx_list,T7_weak_negative_crosstalk_penalty_idx_list]); 
% sigma70_strong_crosstalk_penalty_cutoff = penalty_cutoff([sigma70_strong_positive_crosstalk_penalty_idx_list,sigma70_strong_negative_crosstalk_penalty_idx_list]); 
% sigma70_weak_crosstalk_penalty_cutoff = penalty_cutoff([sigma70_weak_positive_crosstalk_penalty_idx_list,sigma70_weak_negative_crosstalk_penalty_idx_list]); 
% % penalty_cutoff_reorg = [penalty_cutoff(1:6),T7_strong_crosstalk_penalty_cutoff,T7_weak_crosstalk_penalty_cutoff,sigma70_strong_crosstalk_penalty_cutoff,sigma70_weak_crosstalk_penalty_cutoff,...
% %     penalty_cutoff(end-1:end)];
% 
%     % Generate satisfy flag matrix that is #samples * 9
% satisfy_flag_matrix_streamlined = [baseline_dev <= penalty_cutoff(1),crosstalk_ratio_dev <= penalty_cutoff(2),all(baseline_penalty <= penalty_cutoff(3:6),2),...
%     all(T7_strong_crosstalk_penalty <= T7_strong_crosstalk_penalty_cutoff,2),all(T7_weak_crosstalk_penalty <= T7_weak_crosstalk_penalty_cutoff,2),...
%     all(sigma70_strong_crosstalk_penalty <= sigma70_strong_crosstalk_penalty_cutoff,2),all(sigma70_weak_crosstalk_penalty <= sigma70_weak_crosstalk_penalty_cutoff,2),...
%     large_positive_crosstalk_penalty <= penalty_cutoff(end - 1),residual_mRNA_penalty <= penalty_cutoff(end)];
% 
%     % Generate all combinations 
% num_layers = 9; 
% combinations_in_layer = cell(num_layers,1);
% for i = 1:num_layers
%     comb_mat = nchoosek(1:9,i); 
%     combinations_in_layer{i} = comb_mat; 
% 
% end
% 
%     % Calculate number of samples (flux in Sankey diagram)
% satisfied_layer = cell(num_layers,1);
% for i = 1:num_layers
%     comb_mat = combinations_in_layer{i}; 
%     satisfied_layer_node = nan(size(comb_mat,1),1); 
%     for comb_idx = 1:size(comb_mat,1)
%         selected_criteria_idx = comb_mat(comb_idx,:); 
%         num_satisfied = sum(all(satisfy_flag_matrix_streamlined(:,selected_criteria_idx),2)); 
%         satisfied_layer_node(comb_idx) = num_satisfied; 
%     end
%     satisfied_layer{i} = satisfied_layer_node; 
% 
% end

%     % Calculate total number of links
% 
% layers_for_sankey = cell(100000,3);
% layer_add_idx = 1; 
% layer_alphabets = {'A','B','C','D','E','F','G','H','I'}; 
% for layer_idx = 1:num_layers - 1
%     this_layer = combinations_in_layer{layer_idx};
%     next_layer = combinations_in_layer{layer_idx + 1}; 
%     this_layer_letter = layer_alphabets{layer_idx};
%     next_layer_letter = layer_alphabets{layer_idx + 1};
% 
%     for node_idx_this = 1:size(this_layer,1)
%         layer_node_this_idx = this_layer(node_idx_this,:); 
%         for node_idx_next = 1:size(next_layer,1)
%             layer_node_next_idx = next_layer(node_idx_next,:);
%             num_satisfied_layer_node = satisfied_layer{layer_idx + 1}(node_idx_next); 
%             if all(ismember(layer_node_this_idx,layer_node_next_idx))
%                 layers_for_sankey{layer_add_idx,1} = strcat(this_layer_letter,num2str(node_idx_this));
%                 layers_for_sankey{layer_add_idx,2} = strcat(next_layer_letter,num2str(node_idx_next));
%                 layers_for_sankey{layer_add_idx,3} = num_satisfied_layer_node; 
%                 layer_add_idx = layer_add_idx + 1; 
%             end
%         end
%     end
% 
% end
% 
% keep_idx_list = all(~cellfun(@isempty,layers_for_sankey),2); 
% filtered_layers_for_sankey = layers_for_sankey(keep_idx_list,:); 
% 
% figure; 
% SK=SSankey(filtered_layers_for_sankey(:,1),filtered_layers_for_sankey(:,2),filtered_layers_for_sankey(:,3));
% SK.draw()

%     % Calculate number of criteria each row meets
% num_satisfied_criteria_per_sample_streamlined = sum(satisfy_flag_matrix_streamlined,2); 
% num_satisfied_sample_all = nan(size(satisfy_flag_matrix_streamlined,2),1);
% for num_satisfied = 1:size(satisfy_flag_matrix_streamlined,2)
%     num_satisfied_sample = sum(num_satisfied_criteria_per_sample_streamlined >= num_satisfied); 
%     num_satisfied_sample_all(num_satisfied) = num_satisfied_sample; 
% end
% % 
% % labels = {'>= 5','>= 10','>= 15','>= 20','>= 25','>= 30','>= 32'};
% % ,'Labels',labels
% funnelchart( num_satisfied_sample_all)


%% Use criteria above to select sampled parameters for penalty breakdown visualization 
% num_satisfied_8_sample_idx_list = find(num_satisfied_criteria_per_sample_streamlined >= 8); 
% high_res_penalty_terms_oi = sampling_result_file2.all_high_res_penalty_terms(num_satisfied_8_sample_idx_list,:); 
% 
% species_name_list = simFunction_sampling2.Observables.Name;
% species_option_struct = struct('promotor_oi','sigma70_strong','conc_idx',1);
% species_option_struct.species_name_list = species_name_list; 
% 
% % Extract all information related to the sampled parameters of interest 
% idx_oi = num_satisfied_8_sample_idx_list(3);
% sampled_params_oi = init_param_file.satisfied_sampled_params(idx_oi,:); 
% simFunction = simFunction_sampling2; 
% dosing_information = dosing_information_sampling2; 
% problemObject = problemObject_sampling2; 
% save('test_save_files/20241023_model_oi_info.mat','simFunction','dosing_information','sampled_params_oi','problemObject',...
%     'tStart','tEnd'); 

% for k = 3:3%1:length(num_satisfied_8_sample_idx_list)
%     idx_oi = num_satisfied_8_sample_idx_list(k);
%     sampled_params_oi = init_param_file.satisfied_sampled_params(idx_oi,:); 
%     [simulated_time,simulated_data] = simFunction_sampling2(sampled_params_oi,tEnd,dosing_information_sampling2,tStart:tEnd);
% 
%     % plot_simulated(simulated_time,simulated_data,'baseline')
%     % plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
%     plot_simulated(simulated_time,simulated_data,'BundledSpecies',[],[],species_option_struct)
%     plot_simulated(simulated_time,simulated_data,'exhaustive')
% end

% plot_qual_obj_dist(high_res_penalty_terms_oi,penalty_cutoff)
%% Organize sampling results 

% % TestPrototype
% init_param_file = load('test_save_files/20241007_lhs_sampled_params_TestBindingSites.mat');
% all_sampled_params_from_file = init_param_file.satisfied_sampled_params; 
% all_sampled_params = nan(size(all_sampled_params_from_file));
% all_high_res_penalty_terms = nan(size(all_sampled_params_from_file,1),32);
% 
% for iter = 1:385
%     result_file_name = sprintf('param_est_run_save/20241007_TestPrototype_param_sampling_constrained_run%d.mat',iter);
%     if exist(result_file_name,'file')
%         result_file = load(result_file_name);
%         sampled_params = result_file.sampled_params_selected; 
%         high_res_penalty_terms = result_file.high_res_all_penalty_terms; 
%         start_idx = (iter - 1) * 1000 + 1; 
%         end_idx = iter * 1000; 
%         all_sampled_params(start_idx:end_idx,:) = sampled_params; 
% 
%         % Sanity check
%         if ~all(all(sampled_params == all_sampled_params_from_file(start_idx:end_idx,:)))
%             warning('Check sampled parameters')
%         end
% 
%         all_high_res_penalty_terms(start_idx:end_idx,:) = high_res_penalty_terms; 
%     end
% 
% end

% % TestBindingSites
% init_param_file = load('test_save_files/20241007_lhs_sampled_params_TestBindingSites.mat');
% all_sampled_params_from_file = init_param_file.satisfied_sampled_params; 
% problemObject = init_param_file.problemObject;
% dosing_information = init_param_file.dosing_information; 
% all_sampled_params = nan(size(all_sampled_params_from_file));
% all_high_res_penalty_terms = nan(size(all_sampled_params_from_file,1),32);
% 
% for iter = 1:3841
%     result_file_name = sprintf('param_est_run_save/20241007_TestBindingSites_param_sampling_constrained_run%d.mat',iter);
%     if exist(result_file_name,'file')
%         result_file = load(result_file_name);
%         sampled_params = result_file.sampled_params_selected; 
%         high_res_penalty_terms = result_file.high_res_all_penalty_terms; 
%         start_idx = (iter - 1) * 100 + 1; 
%         end_idx = iter * 100; 
%         all_sampled_params(start_idx:end_idx,:) = sampled_params; 
% 
%         % Sanity check
%         if ~all(all(sampled_params == all_sampled_params_from_file(start_idx:end_idx,:)))
%             warning('Check sampled parameters')
%         end
% 
%         all_high_res_penalty_terms(start_idx:end_idx,:) = high_res_penalty_terms; 
%     end
% 
% end
% 
% save('test_save_files/202410_TestBindingSites_sampled.mat','all_sampled_params','all_high_res_penalty_terms','problemObject','dosing_information');

%% Fit txtlsim to show that it cannot capture experimental trends 
%     % Create problemObject
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_run3.xlsx';
% % param_info_path_newMech = 'parameters_v2_newMech.xlsx';
% problemObject = setProblemObject_v3(group_description,[],param_info_path); 
% 
% save_path_id = '_txtlsim_run3';
% wrapper_fit_to_data_w_probObject_v2(problemObject,save_path_id)

% result_file = load('param_est_run_save/20241014_param_est_run1438__txtlsim_run3.mat');
% problemObject = result_file.problemObject; 
% all_fitResults = result_file.all_fitResults; 
% 
% % Organize fitted results 
% allSSE = nan(length(all_fitResults),1);
% all_estimated_params = nan(length(all_fitResults),height(all_fitResults{1}.ParameterEstimates));
% all_simulated_time = cell(length(all_fitResults),1);
% all_simulated_data = cell(length(all_fitResults),1); 
% all_high_res_penalty_term = nan(length(all_fitResults),32); 
% for iter = 1:length(all_fitResults)
%     fitResult = all_fitResults{iter};
%     SSE = fitResult.SSE;
%     estimated_params = fitResult.ParameterEstimates.Estimate; 
%     fitted_data = fitted(fitResult); 
%     [simulated_time,simulated_data] = process_crosstalk_data_from_source(fitted_data,'simulated_data','keep_all');
% 
%     % Assign 
%     allSSE(iter) = SSE; 
%     all_estimated_params(iter,:) = estimated_params; 
%     all_simulated_time{iter} = simulated_time; 
%     all_simulated_data{iter} = simulated_data; 
%     [scaled_penalty_vec,penalty_term_labels,penalty_term_length,~] = wrapper_calculate_obj_higher_resolution_from_data(simulated_time,simulated_data,...
%     problemObject,mode,num_conc,num_promotor);
%     all_high_res_penalty_term(iter,:) = scaled_penalty_vec; 
% 
% end

%% Test visualization of qualitative objectve for parameter sets 


% 
% plot_qual_obj_dist(all_high_res_penalty_term,penalty_cutoff)


%% Test new plot function to bundle bound species
%     % Get a random simulated time and data 
% load('test_save_files/20241007_satisfy_constraints_param_oi.mat');
% 
% species_name_list = simFunction_sampling2.Observables.Name;
% species_option_struct = struct('promotor_oi','sigma70_weak','conc_idx',1);
% species_option_struct.species_name_list = species_name_list; 
% plot_simulated(simulated_time,simulated_data,'BundledSpecies',[],[],species_option_struct)


%% Feasibitliy check for new mechanism 
% missing_idx_list = []; 
% for iter = 1:100
%     result_fileName = sprintf('param_est_run_save/20241010_TestBindingSites_Feasibility_Check_iter%d_updated.mat',iter);
%     if ~exist(result_fileName,'file')
%         wrapper_run_feasibility_check("TestBindingSites",iter); 
%         % missing_idx_list = [missing_idx_list,iter]; 
%     end
% 
% end

%% TX-level crosstalk integration 
    
% 
%     % Load fitted TX-level crosstalk 
% TX_fitted_result_file = load('param_est_run_save/20230509_param_est_run1509_.mat');
% [~,min_idx] = min(TX_fitted_result_file.metric_summary); 
% opt_fitResult = TX_fitted_result_file.all_fitResults{min_idx}; 
% simData = fitted(opt_fitResult);
% [simulated_time,simulated_data] = process_crosstalk_data_from_source(simData,'simulated_data');

    % Set problemObject for TX crosstalk 
% group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject = setProblemObject_v3(group_description,[],param_info_path);
% wrapper_fit_to_data_w_probObject_v2(problemObject,'TX')

%     % Load updated fitted TX-level crosstalk results 
% TX_fitted_result_file_updated = load('param_est_run_save/20241008_param_est_run1356_TX.mat'); 
%     % Compile these fitting results 
% metric_summary_updated = nan(length(TX_fitted_result_file_updated.all_fitResults),1); 
% all_estimated_params = nan(length(TX_fitted_result_file_updated.all_fitResults),31); 
% for iter = 1:length(TX_fitted_result_file_updated.all_fitResults)
% 
%     fitResult = TX_fitted_result_file_updated.all_fitResults{iter}; 
%     estimated_params = [fitResult.ParameterEstimates.Estimate];
%     all_estimated_params(iter,:) = estimated_params'; 
%     metric_summary_updated(iter) = fitResult.SSE; 
% 
%     % % Plot out predicted time-course for every iteration
%     % simData_temp = fitted(fitResult); 
%     % [simulated_time_updated_temp,simulated_data_updated_temp] = ...
%     %     process_crosstalk_data_from_source(simData_temp,'simulated_data');
%     % plot_time_course_crosstalk(simulated_time_updated_temp,simulated_data_updated_temp)
% 
% end

% [~,min_idx_2] = min(metric_summary_updated); 
% opt_fitResult_updated = TX_fitted_result_file_updated.all_fitResults{min_idx_2}; 
% simData_updated = fitted(opt_fitResult_updated);
% [simulated_time_updated,simulated_data_updated] = process_crosstalk_data_from_source(simData_updated,'simulated_data');

% 
% plot_time_course_crosstalk(simulated_time,simulated_data)
% plot_time_course_crosstalk(simulated_time_updated,simulated_data_updated)

%     % Plot out experimental data 
% exp_data = problemObject.Data;
% [exp_time,exp_data] = process_crosstalk_data_from_source(exp_data,'grouped_data');
% plot_time_course_crosstalk(exp_time,exp_data)

    %     % Temp check differences in parameter info
    % parameters_v2_file = readtable('parameters_v2.xlsx');
    % parameters_test2_file = readtable('parameters_test2.xlsx'); 
    % param_name_list_v2 = parameters_v2_file.Name; 
    % param_name_list_test2 = parameters_test2_file.Name; 
    % shared_param_name_list = intersect(param_name_list_test2,param_name_list_v2); 
    % v2_only_param_name_list = setdiff(param_name_list_v2,shared_param_name_list); 
    % test2_only_param_name_list = setdiff(param_name_list_test2,shared_param_name_list); 
    %     % For the shared parameters, check if they have the same initial
    %     % values and LB&UB
    % for share_param_idx = 1:length(shared_param_name_list)
    %     shared_param_name = shared_param_name_list{share_param_idx};
    %     param_idx_v2 = strcmp(param_name_list_v2,shared_param_name);
    %     param_idx_test2 = strcmp(param_name_list_test2,shared_param_name);
    % 
    %     param_vals_v2 = parameters_v2_file{param_idx_v2,["InitVal","LB","UB"]};
    %     param_vals_test2 = parameters_test2_file{param_idx_test2,["InitVal","LB","UB"]};
    %     if ~(param_vals_test2 == param_vals_v2)
    %         fprintf(shared_param_name)
    %         fprintf('\n Test2 values:')
    %         param_vals_test2
    %         fprintf('\n v2 values:')
    %         param_vals_v2
    %     end
    % end

%% Initial Spike Investigation
% load('test_save_files/20241007_initial_spike_crosstalk_ratio_param_oi.mat'); 

% crosstalk_ratio = calculate_crosstalk_ratio_v2(simulated_time,simulated_data,num_conc,num_promotor,mode);

% species_name_list = simFunction_sampling2.Observables.Name;
% species_option_struct_conc1 = struct('promotor_oi','sigma70_strong','conc_idx',1);
% species_option_struct_conc3 = struct('promotor_oi','sigma70_strong','conc_idx',3);
% species_option_struct_conc1.species_name_list = species_name_list; 
% species_option_struct_conc3.species_name_list = species_name_list; 
% plot_simulated(simulated_time,simulated_data,'ConcSpecies',[],[],species_option_struct_conc1)
% plot_simulated(simulated_time,simulated_data,'ConcSpecies',[],[],species_option_struct_conc3)


    % Get information on parameters to be estimated 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% 
% 
% problemObject = setProblemObject_v3(group_description,[],param_info_path); 
% [simulated_time,simulated_data] = process_crosstalk_data_from_source(problemObject.Data,'grouped_data');

% 
% no_empty_conc1_timeVec = simulated_time{1};
% no_empty_conc1_data = simulated_data{1};
% no_empty_conc3_timeVec = simulated_time{3};
% no_empty_conc3_data = simulated_data{3};
% 
% empty_conc1_timeVec = simulated_time{15};
% empty_conc1_data = simulated_data{15};
% empty_conc3_timeVec = simulated_time{17};
% empty_conc3_data = simulated_data{17};
% 
% emptyT7_conc1_timeVec = simulated_time{22};
% emptyT7_conc1_data = simulated_data{22};
% emptyT7_conc3_timeVec = simulated_time{24};
% emptyT7_conc3_data = simulated_data{24};
% 
% figure;
% plot(no_empty_conc1_timeVec,no_empty_conc1_data(:,1),'Color','g','LineWidth',1.5,'LineStyle','-')
% hold on 
% plot(no_empty_conc3_timeVec,no_empty_conc3_data(:,1),'Color','g','LineWidth',1.5,'LineStyle','--')
% plot(empty_conc1_timeVec,empty_conc1_data(:,1),'Color','b','LineWidth',1.5,'LineStyle','-')
% plot(empty_conc3_timeVec,empty_conc3_data(:,1),'Color','b','LineWidth',1.5,'LineStyle','--')
% plot(emptyT7_conc1_timeVec,emptyT7_conc1_data(:,1),'Color','r','LineWidth',1.5,'LineStyle','-')
% plot(emptyT7_conc3_timeVec,emptyT7_conc3_data(:,1),'Color','r','LineWidth',1.5,'LineStyle','--')
% xlabel('Time(s)')
% ylabel('GFP concentration (nM)')
% set(gca,'FontSize',12)
% legend('No Empty 0.5nM','No Empty 2.5nM','Empty sigma70 0.5nM','Empty sigma70 2.5nM','Empty T7 0.5nM','Empty T7 2.5nM')

% Plot out baseline for sigma70 weak lower concentrations 
% figure;
% for data_idx = 85:88
%     simulated_time_single = simulated_time{data_idx};
%     simulated_data_single = simulated_data{data_idx};
% 
%     plot(simulated_time_single,simulated_data_single(:,1),'LineWidth',1.5)
%     hold on 
% 
% end
% xlabel('Time(s)')
% ylabel('GFP concentration (nM)')
% set(gca,'FontSize',12)
% legend('0.5 nM','1 nM','2.5 nM','5 nM')
%% Compile Prototype feasibility check runs results 
% num_iter = 100;
% 
% all_fval = nan(num_iter,1);
% all_init_params = nan(num_iter,31);
% all_est_params = nan(num_iter,31);
% all_exitFlag = nan(num_iter,1);
% all_output = cell(num_iter,1); 
% 
% for iter = 1:num_iter 
% 
%     result_fileName = sprintf('param_est_run_save/202410_TestPrototype_Feasibility_Check_iter%d.mat',iter);
%     if exist(result_fileName,'file')
%         result_file = load(result_fileName); 
%         all_fval(iter) = result_file.fval;
%         all_init_params(iter,:) = result_file.init_params';
%         all_est_params(iter,:) = result_file.opt_params';
%         all_exitFlag(iter) = result_file.exitflag; 
%         all_output{iter} = result_file.output; 
%     else
%         fprintf('Missing iter #%d',iter)
%     end
% 
% end
% 
% save('test_save_files/20241006_TestPrototype_Feasibility_Check_summary_checkpoint1.mat'); 

%% Check fitting vs. sampling results 

    % Load result files 
%         % FITTING (w qual objective function) 
% large_scale_fitting_result_file1 = load('test_save_files/20240529_large_scale_fitting_penalty_term.mat');
% problemObject_fitting = large_scale_fitting_result_file1.problemObject; 
% simFunction_fitting = large_scale_fitting_result_file1.simFunction; 
% dosing_information_fitting = create_dosing_info_from_problemObject(problemObject_fitting); 

        % FITTING 2 (w SSE as obj function) 
% large_scale_fitting_result_file2 = load('param_est_run_save/20231106_large_scale_run.mat'); % Note that this is before the correct toxin mechanism 
% large_scale_fitting_sample_result_file = load('param_est_run_save/20231102_param_est_run1531_1.mat'); 

        %SAMPLING 1 
% sampling_result_file1 = load('test_save_files/202408_high_res_penalty_valid_terms.mat','valid_all_penalty_terms','valid_all_sampled_params','valid_high_res_penalty_terms');
% init_param_file1 = load('test_save_files/20240717_param_init_cond_sampling_constrained.mat','');
% problemObject_sampling1 = init_param_file1.problemObject; 
% all_dosing_information = init_param_file1.all_dosing_information; 
% 
%     % Create simFunction using problemObject
% all_model_species = {problemObject_sampling1.Model.Species.Name}; 
% reporter_species_idx = find(strcmp(all_model_species,reporter_species_name)); 
% additional_track_species = [all_model_species(1:reporter_species_idx - 1),all_model_species(reporter_species_idx + 1:end)];
% simFunction_sampling1 = create_simFun_from_problemObject(problemObject_sampling1,additional_track_species); 

        % %SAMPLING 2
% init_param_file = load('test_save_files/20240819_lhs_sampled_params.mat'); 
% all_satisfied_sampled_params = init_param_file.satisfied_sampled_params; 
% sampling_result_file2 = load('test_save_files/202409_high_res_penalty_sampling_result_summary.mat'); 
% sample_sampling_resultFile2 = load('param_est_run_save/20240819_param_sampling_constrained_run10.mat');
% problemObject_sampling2 = sample_sampling_resultFile2.problemObject; 
% dosing_information_sampling2 = init_param_file.dosing_information; 
%     % Create simFunction using problemObject
% all_model_species2 = {problemObject_sampling2.Model.Species.Name}; 
% reporter_species_idx2 = find(strcmp(all_model_species2,reporter_species_name)); 
% additional_track_species2 = [all_model_species2(1:reporter_species_idx2 - 1),all_model_species2(reporter_species_idx2 + 1:end)];
% simFunction_sampling2 = create_simFun_from_problemObject(problemObject_sampling2,additional_track_species2); 
% 
%     % Use plot_simulated to compare
% num_case = 5; 

%     % FITTING 
% %     % %(use top 5 min obj value)
% % [~,sort_idx] = sort(large_scale_fitting_result_file1.all_obj_value,'ascend');
%     % Use cases with satisfied qualitative phenotypes 
% satisfy_idx_list = (large_scale_fitting_result_file1.all_crosstalk_ratio_transition_penalty < 0) & (large_scale_fitting_result_file1.all_log_baseline_data_penalty < 0);
% satisfied_estimated_params = large_scale_fitting_result_file1.all_estimated_params_mat(satisfy_idx_list,:); 
% satisfied_obj_value = large_scale_fitting_result_file1.all_obj_value(satisfy_idx_list); 
% [~,sort_idx_list] = sort(satisfied_obj_value,'ascend'); 
% 
% for i = 1:num_case
%     selected_idx = sort_idx_list(i);
%     estimated_params_fitting = satisfied_estimated_params(selected_idx,:);
%     [simulated_time,simulated_data] = simFunction_fitting(estimated_params_fitting,tEnd,dosing_information_fitting,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% 
% end
% 
% % i = 5 (selected_idx = 1533) gives a reasonable representation) 

    % FITTING 2
% [~,]

    % SAMPLING 1 
% satisfy_idx_list = sampling_result_file1.valid_high_res_penalty_terms(:,1) < median(log10(sampling_result_file1.valid_all_penalty_terms(:,1)),'omitmissing') &...
%     sampling_result_file1.valid_high_res_penalty_terms(:,2) < median(log10(sampling_result_file1.valid_all_penalty_terms(:,2)),'omitmissing') &...
%     sampling_result_file1.valid_high_res_penalty_terms(:,end) < median(log10(sampling_result_file1.valid_all_penalty_terms(:,end)),'omitmissing') &...
%     sampling_result_file1.valid_all_penalty_terms(:,3) < 0 & sampling_result_file1.valid_all_penalty_terms(:,4) < 0;
% clear init_param_file1
% meta_all_dosing_information = repmat(all_dosing_information,[10000,1]); 
% satisfied_dosing_information = meta_all_dosing_information(satisfy_idx_list);
% satisfied_sampled_params = sampling_result_file1.valid_all_sampled_params(satisfy_idx_list,:); 
% satisfied_high_res_penalty_term = sampling_result_file1.valid_high_res_penalty_terms(satisfy_idx_list,:);
% 
% draft_idx_list = randperm(sum(satisfy_idx_list)); 
% for j = 1:num_case
%     selected_idx = draft_idx_list(j);
%     estimated_params_sampling1 = satisfied_sampled_params(selected_idx,:);
%     dosing_information = meta_all_dosing_information{selected_idx}; 
%     [simulated_time,simulated_data] = simFunction_sampling1(estimated_params_sampling1,tEnd,dosing_information,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% end
    
    % SAMPLING 2   
% satisfy_idx_list = sampling_result_file2.all_high_res_penalty_terms(:,1) < median(log10(sampling_result_file2.all_high_res_penalty_terms(:,1)),'omitmissing') &...
% sampling_result_file2.all_high_res_penalty_terms(:,2) < median(log10(sampling_result_file2.all_high_res_penalty_terms(:,2)),'omitmissing') &...
% sampling_result_file2.all_high_res_penalty_terms(:,end) < median(log10(sampling_result_file2.all_high_res_penalty_terms(:,end)),'omitmissing') &...
% sampling_result_file2.all_high_res_penalty_terms(:,3) < -1e-04 & sampling_result_file2.all_high_res_penalty_terms(:,4) < -1e-04 & ...
% sampling_result_file2.all_high_res_penalty_terms(:,5) < -1e-04 & sampling_result_file2.all_high_res_penalty_terms(:,6) < -1e-04 & ...
% all(sampling_result_file2.all_high_res_penalty_terms(:,19:30) < -0.05,2) & sampling_result_file2.all_high_res_penalty_terms(:,31) < -0.01 &...
% all(sampling_result_file2.all_high_res_penalty_terms(:,7:11) < -0.05,2) & all(sampling_result_file2.all_high_res_penalty_terms(:,13:14) < -0.05,2) & ...
% all(sampling_result_file2.all_high_res_penalty_terms(:,16:18) < -0.05,2);
% satisfy_idx_list = sampling_result_file2.all_high_res_penalty_terms(:,3) < -1e-04 & sampling_result_file2.all_high_res_penalty_terms(:,4) < -1e-04 & ...
% sampling_result_file2.all_high_res_penalty_terms(:,5) < -1e-04 & sampling_result_file2.all_high_res_penalty_terms(:,6) < -1e-04 & ...
% all(sampling_result_file2.all_high_res_penalty_terms(:,19:30) < -0.05,2) & sampling_result_file2.all_high_res_penalty_terms(:,31) < -0.01 &...
% all(sampling_result_file2.all_high_res_penalty_terms(:,7:11) < -0.05,2) & all(sampling_result_file2.all_high_res_penalty_terms(:,13:14) < -0.05,2) & ...
% all(sampling_result_file2.all_high_res_penalty_terms(:,16:18) < -0.05,2);
% satisfied_sampled_params = init_param_file.satisfied_sampled_params(satisfy_idx_list,:); 
% satisfied_high_res_penalty_term = sampling_result_file2.all_high_res_penalty_terms(satisfy_idx_list,:);
% 
% draft_idx_list = randperm(sum(satisfy_idx_list)); 
% for j = 1:6% num_case
%     selected_idx = draft_idx_list(j);
%     estimated_params_sampling2 = satisfied_sampled_params(selected_idx,:);
%     [simulated_time,simulated_data] = simFunction_sampling2(estimated_params_sampling2,tEnd,dosing_information_sampling2,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% end

% Find cases with both baseline deviation and crosstalk ratio deviation
% ranked top
% [~,sort_idx_baseline_dev] = sort(sampling_result_file2.all_high_res_penalty_terms(:,1),'ascend');
% [~,sort_idx_crosstalk_ratio_dev] = sort(sampling_result_file2.all_high_res_penalty_terms(:,2),'ascend');
% [~,rank_baseline_dev] = sort(sort_idx_baseline_dev);
% [~,rank_crosstalk_ratio_dev] = sort(sort_idx_crosstalk_ratio_dev); 
% idx_list_oi = intersect(rank_baseline_dev(1:1000),rank_crosstalk_ratio_dev(1:1000));
% 
% for j = 1:length(idx_list_oi)
%     selected_idx = idx_list_oi(j);
%     estimated_params_sampling2 = all_satisfied_sampled_params(selected_idx,:);
%     [simulated_time,simulated_data] = simFunction_sampling2(estimated_params_sampling2,tEnd,dosing_information_sampling2,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% end


% % Find cases where there is an initial spike 
% idx_oi = draft_idx_list(3); 
% sampled_params_oi = satisfied_sampled_params(3,:);
% [simulated_time,simulated_data] = simFunction_sampling2(sampled_params_oi,tEnd,dosing_information_sampling2,tStart:tEnd);
% save('test_save_files/20241007_initial_spike_crosstalk_ratio_param_oi.mat','dosing_information_sampling2','problemObject_sampling2',...
%     'sampled_params_oi','simFunction_sampling2','simulated_data','simulated_time'); 

% sampled_params_oi = satisfied_sampled_params(draft_idx_list(4),:); 
% [simulated_time,simulated_data] = simFunction_sampling2(sampled_params_oi,tEnd,dosing_information_sampling2,tStart:tEnd);
% plot_simulated(simulated_time,simulated_data,'baseline')
% plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% save('test_save_files/20241007_satisfy_constraints_param_oi.mat','sampled_params_oi','simFunction_sampling2','problemObject_sampling2','dosing_information_sampling2',...
%     "simulated_data",'simulated_time'); 
% Investigae this parameter sets, plot out all species 
    % Checkout T7 strong lowest concentration 
% load('test_save_files/20241007_satisfy_constraints_param_oi.mat');
% 
% species_name_list = simFunction_sampling2.Observables.Name;
% species_option_struct = struct('promotor_oi','sigma70_weak','conc_idx',1);
% species_option_struct.species_name_list = species_name_list; 
% plot_simulated(simulated_time,simulated_data,'ConcSpecies',[],[],species_option_struct)


%% Visual inspection of broc & protein data 
% data_structure_AP_file = load('Data/data_structures/simbio_data_table_updated_AP.mat'); 
% data_structure_FP_file = load('Data/data_structures/simbio_data_table_updated_FP.mat'); 
% data_structure_PE_FP_file = load('Data/data_structures/simbio_data_table_PE_updated_FP.mat'); 


%% For group meeting: Find the current best case for simulation 
    % Also need to do a sanity check on whether the sampled parameters give
    % the correct time-course calculation 
% 
% sampling_result_file = load('test_save_files/202408_high_res_penalty_valid_terms.mat','valid_high_res_penalty_terms','valid_all_sampled_params'); 
% high_res_penalty_terms = sampling_result_file.valid_high_res_penalty_terms; 
% all_sampled_params = sampling_result_file.valid_all_sampled_params;
% 
% sample_result_file = load('param_est_run_save/20240819_param_sampling_constrained_run1.mat'); 
% 
% 
% sample_qual_obj_result_file = load('param_est_run_save/20240822_lhs_sampled_params_opt_run1.mat'); 
% simFunction = sample_qual_obj_result_file.simFunction; 
% dosing_information = create_dosing_info_from_problemObject(sample_qual_obj_result_file.problemObject);
% 
% median_baseline_penalty_value = median(high_res_penalty_terms(:,1));
% median_crosstalk_ratio_penalty_value = median(high_res_penalty_terms(:,2));
% median_other_penalty_value = median(high_res_penalty_terms(:,end)); 
% 
% positive_crosstalk_goals = zeros(1,12); % -0.01 .* ones(1,12);
% negative_crosstalk_goals = zeros(1,12); % -0.01 .* ones(1,12); 
% penalty_goals = [median_baseline_penalty_value,median_crosstalk_ratio_penalty_value,...
%     0,0,0,0,positive_crosstalk_goals,negative_crosstalk_goals,100,...
%     median_other_penalty_value];
% 
% satisfy_high_res_mat = high_res_penalty_terms < penalty_goals; 
% satisfy_idx_list = find(all(satisfy_high_res_mat,2)); 
% 
%     % Plot out simulation for 5 cases
% for idx = 1:1%length(satisfy_idx_list)
%     satisfy_idx = satisfy_idx_list(idx);
%     estimated_params = all_sampled_params(satisfy_idx,:);
% 
%     [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% end

%% Test New Mechanism for Entire Model 
%     % Create problemObject
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path_ctrl = 'parameters_v2.xlsx';
% param_info_path_newMech = 'parameters_v2_newMech.xlsx';
% problemObject_ctrl = setProblemObject_v3(group_description,[],param_info_path_ctrl); 
% problemObject_BindingSite = setProblemObject_v3_TestBindingSite(group_description,[],param_info_path_newMech); 
% problemObject_Prototype = setProblemObject_v3_TestPrototype(group_description,[],param_info_path_newMech); 
% 
% estimated_param_name_list_ctrl = {problemObject_ctrl.Estimated.Name};
% estimated_param_name_list_BindingSite = {problemObject_BindingSite.Estimated.Name};
% estimated_param_name_list_Prototype = {problemObject_Prototype.Estimated.Name};
% 
    % Perform a test run with newly added mRNA degradation mechanism 
% wrapper_fit_to_data_w_probObject(problemObject_ctrl,'ctrl');
% wrapper_fit_to_data_w_probObject(problemObject_BindingSite,'TestBindSites');
% wrapper_fit_to_data_w_probObject(problemObject_Prototype,'TestPrototype');
% 
    % Load result from test runs 
% result_file_ctrl = load('param_est_run_save/20230926_param_est_run1324_ctrl.mat'); 
% result_file_TestPrototype = load('param_est_run_save/20231002_param_est_run1653_TestPrototype.mat'); 
% result_file_TestBindingSite = load('param_est_run_save/20231001_param_est_run1027_TestBindSites.mat'); 
% 
% fitResult_ctrl = result_file_ctrl.all_fitResults{1,1};
% fitResult_TestPrototype = result_file_TestPrototype.all_fitResults{1,1};
% fitResult_TestBindingSite = result_file_TestBindingSite.all_fitResults{1,1}; 
% 
% fitted_data_ctrl = fitted(fitResult_ctrl);
% fitted_data_TestPrototype = fitted(fitResult_TestPrototype);
% fitted_data_TestBindingSite = fitted(fitResult_TestBindingSite); 
% 
% [simulated_time_ctrl,simulated_data_ctrl] = process_crosstalk_data_from_source(fitted_data_ctrl,'simulated_data');
% [simulated_time_TestPrototype,simulated_data_TestPrototype] = process_crosstalk_data_from_source(fitted_data_TestPrototype,'simulated_data');
% [simulated_time_TestBindingSite,simulated_data_TestBindingSite] = process_crosstalk_data_from_source(fitted_data_TestBindingSite,'simulated_data');
% 
% plot_simulated(simulated_time_ctrl,simulated_data_ctrl,'baseline')
% plot_simulated(simulated_time_ctrl,simulated_data_ctrl,'crosstalk_ratio')
% 
% plot_simulated(simulated_time_TestPrototype,simulated_data_TestPrototype,'baseline')
% plot_simulated(simulated_time_TestPrototype,simulated_data_TestPrototype,'crosstalk_ratio')
% 
% plot_simulated(simulated_time_TestBindingSite,simulated_data_TestBindingSite,'baseline')
% plot_simulated(simulated_time_TestBindingSite,simulated_data_TestBindingSite,'crosstalk_ratio')
% 
% figure;
% simData_single_lowConc_no_empty = fitted_data_TestPrototype(1);
% simData_single_highConc_no_empty = fitted_data_TestPrototype(7); 
% simData_single_lowConc_empty = fitted_data_TestPrototype(8);
% simData_single_highConc_empty = fitted_data_TestPrototype(14);
% 
% for species_idx = 76:length(simData_single_lowConc_no_empty.DataNames)
% 
%     subplot(5,5,species_idx - 75)
% 
%     plot(simData_single_highConc_no_empty.Time,simData_single_highConc_no_empty.Data(:,species_idx),'LineWidth',1.5,'Color','g')
%     hold on 
%     plot(simData_single_highConc_empty.Time,simData_single_highConc_empty.Data(:,species_idx),'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%     title(simData_single_highConc_no_empty.DataNames{species_idx}); 
% 
% end
% 
% % Blank test space
% test_toy_model = sbiomodel('simpleModel');
% test_binding_rxn = addreaction(test_toy_model,'A + B <-> C');
% kineticLaw = addkineticlaw(test_binding_rxn,'MassAction');
% binding_F = addparameter(kineticLaw,'k_f',0.5);
% binding_R = addparameter(kineticLaw,'k_r',0.5);
% kineticLaw.ParameterVariableNames = {'k_f','k_r'};
% 
%     % In this test toy model, I want to create a new reaction and its
%     % kineticLaw.Parameters contain binding_F and binding_R without
%     % assuming that binding_F & binding_R exist in working space 
% binding_parameters = get(test_binding_rxn.KineticLaw,'Parameters'); 
% binding_F_2 = binding_parameters(1); 
% binding_R_2 = binding_parameters(2); 
% 
% test_new_rxn = addreaction(test_toy_model,'D + E <-> F');
% kineticLaw_new_rxn = addkineticlaw(test_new_rxn,'MassAction');
% binding_F_2_for_rxn = addparameter(kineticLaw_new_rxn,binding_F_2.Name,binding_F_2.Value);
% binding_R_2_for_rxn = addparameter(kineticLaw_new_rxn,binding_R_2.Name,binding_R_2.Value);
% set(kineticLaw_new_rxn, 'ParameterVariableNames', {binding_parameters.Name});



%% Analyze optimization runs with qual obj  
% load('test_save_files/20240919_lhd_sampled_params_optimizatio_qual_obj_result_summary.mat');

% % Plot out penalty distribution
% figure; 
% for penalty_idx = 1:size(all_final_penalty,2)
% 
%     subplot(ceil(sqrt(size(all_final_penalty,2))),ceil(sqrt(size(all_final_penalty,2))),penalty_idx)
% 
%     penalty_label = penalty_term_label{penalty_idx};  
%     init_penalty_term = all_init_penalty(:,penalty_idx);
%     final_penalty_term = all_final_penalty(:,penalty_idx); 
% 
%     % Remove outliers 
%     filtered_init_penalty_term = rmoutliers(init_penalty_term); 
%     filtered_final_penalty_term = rmoutliers(final_penalty_term); 
% 
%     if penalty_idx <= 2 || isequal(penalty_idx,size(all_final_penalty,2))
%         init_penalty_values_for_plot = log10(filtered_init_penalty_term); 
%         final_penalty_values_for_plot = log10(filtered_final_penalty_term); 
%         title_name = strcat('log10',strrep(penalty_label,'_',' ')); 
%     else
%         init_penalty_values_for_plot = filtered_init_penalty_term; 
%         final_penalty_values_for_plot = filtered_final_penalty_term; 
%         title_name = strrep(penalty_label,'_',' '); 
%     end
% 
%     histogram(init_penalty_values_for_plot)
%     hold on 
%     histogram(final_penalty_values_for_plot)
%     ylabel('# Occurence')
%     title(title_name); 
%     set(gca,'FontSize',12);
% 
% end
% legend('Initial Penalty','Final Penalty')

% % Deeper dive into runs that satisfy all constraints 
% sampling_result_file = load('test_save_files/202408_high_res_penalty_valid_terms.mat','valid_high_res_penalty_terms'); 
% high_res_penalty_terms = sampling_result_file.valid_high_res_penalty_terms; 
% 
% median_baseline_penalty_value = median(high_res_penalty_terms(:,1));
% median_crosstalk_ratio_penalty_value = median(high_res_penalty_terms(:,2));
% median_other_penalty_value = median(high_res_penalty_terms(:,end)); 
% 
% % Is there a case where the estimated parameters satisfy all the
% % constraints? 
% satisfy_idx_list = find((log10(all_final_penalty(:,1)) < median_baseline_penalty_value) & (log10(all_final_penalty(:,2)) < median_crosstalk_ratio_penalty_value) ...
%     & all_final_penalty(:,3) < 0 & all_final_penalty(:,4) < -0.1 & all_final_penalty(:,5) < -0.1 & (log10(all_final_penalty(:,7)) < median_other_penalty_value));
% 
% % For each of these cases, check the simulated data 
%     % Get simFunction
% sample_qual_obj_result_file = load('param_est_run_save/20240822_lhs_sampled_params_opt_run1.mat'); 
% simFunction = sample_qual_obj_result_file.simFunction; 
% dosing_information = create_dosing_info_from_problemObject(sample_qual_obj_result_file.problemObject);
% 
%     % Plot out simulation for 5 cases
% for idx = 6:19
%     satisfy_idx = satisfy_idx_list(idx);
%     estimated_params = all_estimated_params(satisfy_idx,:);
% 
%     [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% end

%% Analyze feasibility results 
% feasibilty_check_result_file_run1 = load('test_save_files/20240916_fmincon_feasibility_test_run1.mat'); 
% feasibilty_check_result_file_run2 = load('test_save_files/20240916_fmincon_feasibility_test_run2_feasibilityModeOn.mat');
% feasibilty_check_result_file_run3 = load('test_save_files/20240918_fmincon_feasibility_test_run3_feasibilityModeOn_MultipleStart.mat');
% feasibility_check_result_file_run5 = load('test_save_files/20240923_fmincon_feasibility_test_run5_feasibilityModeOn_MultipleStart.mat'); 


%% Analyze single qualitative objective function 

% qualitative_obj_string_file = load('qualitative_trend_obj_string.mat'); 
% all_obj_string = qualitative_obj_string_file.all_obj_string;
% 
% high_res_example_result_file = load('param_est_run_save/20240821_param_sampling_constrained_run1.mat');
% penalty_term_labels = high_res_example_result_file.penalty_term_labels;
% penalty_term_length = high_res_example_result_file.penalty_term_length; 
% 
% sampling_result_file = load('test_save_files/202408_high_res_penalty_valid_terms.mat','valid_high_res_penalty_terms'); 
% high_res_penalty_terms = sampling_result_file.valid_high_res_penalty_terms; 
% 
% median_baseline_penalty_value = median(high_res_penalty_terms(:,1));
% median_crosstalk_ratio_penalty_value = median(high_res_penalty_terms(:,2));
% median_other_penalty_value = median(high_res_penalty_terms(:,end)); 
% 
% positive_crosstalk_goals = -0.05 .* ones(1,12);
% negative_crosstalk_goals = -0.05 .* ones(1,12); 
% penalty_goals = [median_baseline_penalty_value,median_crosstalk_ratio_penalty_value,-1e-04,-1e-04,-1e-04,-1e-04,positive_crosstalk_goals,negative_crosstalk_goals,-1e-04,median_other_penalty_value];
% 
% 
% figure; 
% for obj_idx = 1:32
% 
%     result_fileName = sprintf('param_est_run_save/20240823_fit_single_qual_obj%d_summary.mat',obj_idx); 
%     result_file = load(result_fileName); 
% 
%     % Extract info from result file struct 
%     all_fval = result_file.all_fval; 
%     all_estimated_params = result_file.all_estimated_params;
%     all_exitFlag = result_file.all_exitFlag; 
%     all_output = result_file.all_output; 
%     problemObject = result_file.problemObject;
%     dosing_information = result_file.dosing_information; 
% 
%     % Plot out fval distribution 
%     subplot(6,6,obj_idx)
%         % Get penalty label for plot title
%     if obj_idx <= 6 
%         penalty_label = penalty_term_labels{obj_idx};
%     elseif obj_idx <= 18
%         penalty_label = strcat(penalty_term_labels{7},num2str(obj_idx - 6));
%     elseif obj_idx <= 30
%         penalty_label = strcat(penalty_term_labels{8},num2str(obj_idx - 18));
%     else
%         penalty_label = penalty_term_labels{obj_idx - 22};  
%     end
%     filtered_all_fval = rmoutliers(all_fval); 
%     histogram(filtered_all_fval)
%     hold on
%     xline(penalty_goals(obj_idx),'LineWidth',1.5,'Color','r'); 
%     title(strrep(penalty_label,'_',' '))
% 
% 
% end

%     % Single objective 
% 
%     % Load a sample result file 
% sample_result_file = load('param_est_run_save/20240823_fit_single_qual_obj1_iter1.mat');
% problemObject = sample_result_file.problemObject;
% dosing_information = sample_result_file.dosing_information; 
% num_params = length(problemObject.Estimated); 
% 
% num_iter = 48; 
% 
% qualitative_obj_string_file = load('qualitative_trend_obj_string.mat'); 
% all_obj_string = qualitative_obj_string_file.all_obj_string;
% for obj_idx = 1:32
% 
%     all_estimated_params = nan(num_iter,num_params); 
%     all_fval = nan(num_iter,1); 
%     all_exitFlag = nan(num_iter,1);
%     all_output = cell(num_iter,1); 
% 
%     for iter = 1:num_iter
%         try 
%             single_obj_opt_result_file_name = sprintf('param_est_run_save/20240823_fit_single_qual_obj%d_iter%d.mat',obj_idx,iter);
%             if ~exist(single_obj_opt_result_file_name,'file')
%                 fprintf('Running obj #32 iter #%d',iter)
%                 run_fit_qualitative_characteristics_single_v2(32,iter)
%             end
%             single_obj_opt_result_file = load(single_obj_opt_result_file_name); 
% 
%             all_estimated_params(iter,:) = single_obj_opt_result_file.estimated_params;
%             all_fval(iter) = single_obj_opt_result_file.fval;
%             all_exitFlag(iter) = single_obj_opt_result_file.exitflag;
%             all_output{iter} = single_obj_opt_result_file.output;
%         catch
%             fprintf('Missing obj %d iter %d',obj_idx,iter)
%         end
% 
%     end
% 
%     save(sprintf('param_est_run_save/20240823_fit_single_qual_obj%d_summary.mat',obj_idx),'problemObject','all_output','all_exitFlag','all_estimated_params',...
%         'all_fval','problemObject','dosing_information')
% end

%% Analyze fitted RNA deg data 
% result_file_imputed = load('test_save_files/20240912_alterantive_kinetics_fitting_imputed_data_expanded_bound.mat'); 
% result_file_original = load('test_save_files/20240912_alterantive_kinetics_fitting_original_data_expanded_bound.mat'); 
% result_file_simulated = load('test_save_files/20240912_alterantive_kinetics_fitting_simulated_data_expanded_bound.mat'); 

%% Compile large-scale sampling results with expanded parameter bounds 
%     % Test load result file 
% result_fileName_example = 'param_est_run_save/20240821_param_sampling_constrained_run856.mat';
% result_file_example = load(result_fileName_example);
% 
% init_param_file = load('test_save_files/20240819_lhs_sampled_params.mat'); 
% num_result_file = ceil(size(init_param_file.satisfied_sampled_params,1) / 100); 
% all_sampled_params = init_param_file.satisfied_sampled_params; 
% 
% 
% % Preassign space for penalty terms 
% all_high_res_penalty_term = nan(size(init_param_file.satisfied_sampled_params,1),size(result_file_example.high_res_all_penalty_terms,2));
% for iter = 1:num_result_file
% 
%     % Load result file
%     result_file_name = sprintf('param_est_run_save/20240821_param_sampling_constrained_run%d.mat',iter); 
%     result_file = load(result_file_name); 
% 
%     % Extract information 
%     sampled_params_selected = result_file.sampled_params_selected;
%     sampled_penalty = result_file.high_res_all_penalty_terms; 
%     if ~exist('problemObject','var')
%         problemObject = result_file.problemObject;
%         dosing_information = result_file.dosing_information;
%         penalty_term_labels = result_file.penalty_term_labels;
%         penalty_term_length = result_file.penalty_term_length; 
%     end
% 
%     % Sanity check for sampled parameters alignment 
%     init_param_file_sampled_params = init_param_file.satisfied_sampled_params((iter - 1) * 100 + 1:iter * 100,:); 
%     if ~all(init_param_file_sampled_params == sampled_params_selected)
%         error('Mismatched sampled parameters')
%     end
% 
%     % Compile 
%     all_high_res_penalty_term((iter - 1) * 100 + 1:iter * 100,:) = sampled_penalty; 
% 
% 
% end

%% Initial analysis of large-scale fitting & single-objective fitting 
    % Large-scale fitting 
% 
%         % SSE
% fitting_result_file_name_SSE_sample = 'param_est_run_save/20240822_lhs_sampled_params_opt_run2.mat'; 
% fitting_result_file_SSE_sample = load(fitting_result_file_name_SSE_sample); 
% 
%         % Qualitative objective 
% fitting_result_file_name_qual_obj_sample = 'param_est_run_save/20240822_lhs_sampled_params_opt_qual_obj_run2.mat'; 
% fitting_result_file_qual_obj_sample = load(fitting_result_file_name_qual_obj_sample); 
% 
% init_params_file = load('test_save_files/20240822_lhs_sampled_params_optimization.mat');
% init_params_file_all_sampled_params = init_params_file.satisfied_sampled_params; 
% 
%         % Check how many iterations are run in each case 
% num_completed_run_SSE = nan(16,1);
% num_completed_run_qual_obj = nan(16,1); 
% for idx = 1:16
%     fitting_result_file_name_SSE = sprintf('param_est_run_save/20240822_lhs_sampled_params_opt_run%d.mat',idx); 
%     fitting_result_file_name_qual_obj = sprintf('param_est_run_save/20240822_lhs_sampled_params_opt_qual_obj_run%d.mat',idx);
%     fitting_result_file_SSE = load(fitting_result_file_name_SSE);
%     fitting_result_file_qual_obj = load(fitting_result_file_name_qual_obj); 
% 
%     num_completed_run_qual_obj(idx) = nnz(~cellfun(@isempty,fitting_result_file_qual_obj.all_fitResults));
%     num_completed_run_SSE(idx) = nnz(~cellfun(@isempty,fitting_result_file_SSE.all_fitResults));
% 
% end
% 
%     % Compile results for results using qualitative features as objective function
%         % Record: init & final params + penalty, fval, 
% num_params = length(fitting_result_file_SSE_sample.all_fitResults{1,1}.estimated_params); 
% all_init_params = nan(16 * 100,num_params);
% all_estimated_params = nan(16 * 100,num_params);
% all_optimal_obj = nan(16 * 100,1); 
% all_init_penalty = nan(16 * 100,7);
% all_final_penalty = nan(16 * 100,7); 
% all_exitFlag = nan(16 * 100,1); 
% 
% for idx = 1:16
%         % Saving names were mixed up - this is counterintuitive but correct
%     % fitting_result_file_name_qual_obj = sprintf('param_est_run_save/20240822_lhs_sampled_params_opt_run%d.mat',idx); 
%     fitting_result_file_name_qual_obj = sprintf('param_est_run_save/20240822_lhs_sampled_params_opt_qual_obj_run%d.mat',idx);
%     fitting_result_file_qual_obj= load(fitting_result_file_name_qual_obj);
%     % fitting_result_file_qual_obj = load(fitting_result_file_name_qual_obj); 
% 
%     % % Sanity check with initial parameter file 
%     % sampled_params_from_init_param_file = init_params_file_all_sampled_params((idx - 1) * 100 + 1:idx * 100,:);
%     % sampled_params_from_result_file = fitting_result_file_qual_obj.; 
%     % if ~all(all(sampled_params_from_init_param_file == sampled_params_from_result_file))
%     %     error('Sanity check on initial parameter failed')
%     % end
% 
%     % Access fields in result file 
%     all_fitResults = fitting_result_file_qual_obj.all_fitResults;
%     estimated_params_org = nan(length(all_fitResults),num_params);
%     fval_org = nan(length(all_fitResults),1); 
%     exitflag_org = nan(length(all_fitResults),1); 
% 
%     for iter = 1:length(all_fitResults)
%         estimated_params_org(iter,:) = all_fitResults{iter}.estimated_params; 
%         fval_org(iter) = all_fitResults{iter}.fval;
%         exitflag_org(iter) = all_fitResults{iter}.exitflag; 
%     end
% 
%     % Assign values 
%     all_init_params((idx - 1) * 100 + 1:idx * 100,:) = init_params_file_all_sampled_params((idx - 1) * 100 + 1:idx * 100,:); 
%     all_estimated_params((idx - 1) * 100 + 1:idx * 100,:) = estimated_params_org;
%     all_optimal_obj((idx - 1) * 100 + 1:idx * 100,:) = fval_org; 
%     all_init_penalty((idx - 1) * 100 + 1:idx * 100,:) = fitting_result_file_qual_obj.all_penalty_term_init; 
%     all_final_penalty((idx - 1) * 100 + 1:idx * 100,:) = fitting_result_file_qual_obj.all_penalty_term_final; 
%     all_exitFlag((idx - 1) * 100 + 1:idx * 100,:) = exitflag_org;
% 
% end
% save('test_save_files/20240919_lhd_sampled_params_optimizatio_SSE_result_summary.mat','all_init_params','all_estimated_params','all_init_penalty',...
%     'all_final_penalty','all_exitFlag');

%     % Single objective 
% 
%     % Load a sample result file 
% sample_result_file = load('param_est_run_save/20240823_fit_single_qual_obj1_iter1.mat');
% problemObject = sample_result_file.problemObject;
% dosing_information = sample_result_file.dosing_information; 
% num_params = length(problemObject.Estimated); 
% 
% num_iter = 48; 
% 
% qualitative_obj_string_file = load('qualitative_trend_obj_string.mat'); 
% all_obj_string = qualitative_obj_string_file.all_obj_string;
% for obj_idx = 1:9
% 
%     all_estimated_params = nan(num_iter,num_params); 
%     all_fval = nan(num_iter,1); 
%     all_exitFlag = nan(num_iter,1);
%     all_output = cell(num_iter,1); 
% 
%     for iter = 1:num_iter
%         try 
%             single_obj_opt_result_file_name = sprintf('param_est_run_save/20240823_fit_single_qual_obj%d_iter%d.mat',obj_idx,iter);
%             single_obj_opt_result_file = load(single_obj_opt_result_file_name); 
% 
%             all_estimated_params(iter,:) = single_obj_opt_result_file.estimated_params;
%             all_fval(iter) = single_obj_opt_result_file.fval;
%             all_exitFlag(iter) = single_obj_opt_result_file.exitflag;
%             all_output{iter} = single_obj_opt_result_file.output;
%         catch
%             fprintf('Missing obj %d iter %d',obj_idx,iter)
%         end
% 
%     end
% 
%     save(sprintf('param_est_run_save/20240823_fit_single_qual_obj%d_summary.mat',obj_idx),'problemObject','all_output','all_exitFlag','all_estimated_params',...
%         'all_fval','problemObject','dosing_information')
% end
% 
%     % Weights 
% result_file_name_w_weights = 'param_est_run_save/20240829_param_est_run1207_w_weight.mat';
% result_file_name_wo_weights = 'param_est_run_save/20240829_param_est_run1207_wo_weight.mat';
% result_file_w_weights = load(result_file_name_w_weights);
% result_file_wo_weights = load(result_file_name_wo_weights); 
% 
% % Let's analyze the weights results with whtat we have so far
% num_completed_run_w_weights = nnz(~cellfun(@isempty,result_file_w_weights.all_fitResults));
% num_completed_run_wo_weights = nnz(~cellfun(@isempty,result_file_wo_weights.all_fitResults));
% 
%     % First column ctrl (wo weights), 2nd column w weights 
%      % Also need to check whether the phenotypes are captured - calculate high resolution penalty values 
% SSE_comparison = nan(min([num_completed_run_wo_weights,num_completed_run_w_weights]),2); 
% all_high_res_penalty_w_weights = nan(min([num_completed_run_wo_weights,num_completed_run_w_weights]),32); 
% all_high_res_penalty_wo_weights = nan(min([num_completed_run_wo_weights,num_completed_run_w_weights]),32); 
% for iter = 55:min([num_completed_run_wo_weights,num_completed_run_w_weights]) 
% 
%     if isequal(iter,1)
%         problemObject = result_file_w_weights.problemObject; 
%         dosing_information = create_dosing_info_from_problemObject(problemObject); 
%     end
% 
%     fitResult_w_weights = result_file_w_weights.all_fitResults{iter};
%     fitResult_wo_weights = result_file_wo_weights.all_fitResults{iter};
% 
%     estimated_params_w_weights = [fitResult_w_weights.ParameterEstimates.Estimate]';
%     estimated_params_wo_weights = [fitResult_wo_weights.ParameterEstimates.Estimate]';
% 
%     SSE_comparison(iter,:) = [fitResult_wo_weights.SSE,fitResult_w_weights.SSE];
% 
%     [scaled_penalty_vec_w_weights,penalty_term_labels,~,~] = wrapper_calculate_obj_higher_resolution(estimated_params_w_weights,...
%         problemObject,dosing_information,[],[],mode,num_conc,num_promotor);
%     [scaled_penalty_vec_wo_weights,~,~,~] = wrapper_calculate_obj_higher_resolution(estimated_params_wo_weights,...
%         problemObject,dosing_information,[],[],mode,num_conc,num_promotor);
% 
%     all_high_res_penalty_w_weights(iter,:) = scaled_penalty_vec_w_weights; 
%     all_high_res_penalty_wo_weights(iter,:) = scaled_penalty_vec_wo_weights;
% 
% end
% 
% % Plot out comparison of high resolution penalty
% for penalty_idx = 1:size(all_high_res_penalty_w_weights,2)
% 
%     subplot(ceil(sqrt(size(all_high_res_penalty_w_weights,2))),ceil(sqrt(size(all_high_res_penalty_w_weights,2))),penalty_idx)
%     bar([all_high_res_penalty_w_weights(:,penalty_idx),all_high_res_penalty_wo_weights(:,penalty_idx)]);
% 
%     if penalty_idx <= 6 
%         penalty_label = penalty_term_labels{penalty_idx};
%     elseif penalty_idx <= 18
%         penalty_label = strcat(penalty_term_labels{7},num2str(penalty_idx - 6));
%     elseif penalty_idx <= 30
%         penalty_label = strcat(penalty_term_labels{8},num2str(penalty_idx - 18));
%     else
%         penalty_label = penalty_term_labels{penalty_idx - 22};  
%     end
%     title(strrep(penalty_label,'_',' '));
% end
% 
% legend('Ctrl','With normalized data')
% 
% % Check how many cases of positive -> negative crosstalk ratio transition
% % were captured
% 
%     % Check whether toxin mechanism is correct 
% num_toxin_captured_w_weights = sum(all_high_res_penalty_w_weights(:,3) < 0 & all_high_res_penalty_w_weights(:,4) < 0 ...
%     & all_high_res_penalty_w_weights(:,5) < 0 & all_high_res_penalty_w_weights(:,6) < 0);
% num_toxin_captured_wo_weights = sum(all_high_res_penalty_wo_weights(:,3) < 0 & all_high_res_penalty_wo_weights(:,4) < 0 ...
%     & all_high_res_penalty_wo_weights(:,5) < 0 & all_high_res_penalty_wo_weights(:,6) < 0);
% num_toxin_decrease_captured_w_weights = sum(all_high_res_penalty_w_weights(:,3) < 0 & all_high_res_penalty_w_weights(:,4) < 0 ...
%     & all_high_res_penalty_w_weights(:,5) < 0);
% num_toxin_decrease_captured_wo_weights = sum(all_high_res_penalty_wo_weights(:,3) < 0 & all_high_res_penalty_wo_weights(:,4) < 0 ...
%     & all_high_res_penalty_wo_weights(:,5) < 0);
% 
%     % Check whether positive -> negative crosstalk ratio transition was
%     % captured 
% num_crosstalk_transition_captured_w_weights = sum(all(all_high_res_penalty_w_weights(:,7:30) < 0,2));
% num_crosstalk_transition_captured_wo_weights = sum(all(all_high_res_penalty_wo_weights(:,7:30) < 0,2));
% num_positive_crosstalk_captured_w_weights = sum(all(all_high_res_penalty_w_weights(:,7:18) < 0,2));
% num_positive_crosstalk_captured_wo_weights = sum(all(all_high_res_penalty_wo_weights(:,7:18) < 0,2));
% num_negative_crosstalk_captured_w_weights = sum(all(all_high_res_penalty_w_weights(:,19:30) < 0,2));
% num_negative_crosstalk_captured_wo_weights = sum(all(all_high_res_penalty_wo_weights(:,19:30) < 0,2));
% 
% figure; 
% xlabels = {'# Toxin Mechanism Captured','# Decrease Toxin Mechanism Captured','# Crosstalk Transition Captured',...
%     '# Positive Crosstalk Captured','# Negative Crosstalk Captured'}; 
% xlabels_for_plot = categorical(xlabels);
% xlabels_for_plot = reordercats(xlabels_for_plot,xlabels);
% bar(xlabels_for_plot,[num_toxin_captured_wo_weights,num_toxin_captured_w_weights;num_toxin_decrease_captured_wo_weights,num_toxin_decrease_captured_w_weights;...
%     num_crosstalk_transition_captured_wo_weights,num_crosstalk_transition_captured_w_weights;num_positive_crosstalk_captured_wo_weights,num_positive_crosstalk_captured_w_weights;...
%     num_negative_crosstalk_captured_wo_weights,num_negative_crosstalk_captured_w_weights]);
% legend('Ctrl','With normalized data')

%% OAT sensitivity analysis for rational parameter bounds relaxation 

% % Bound relaxation for sensitive parameters only 
% 
%     % Load LHS-PRCC results 
% result_file = load('test_save_files/202408_high_res_penalty_valid_terms.mat','valid_high_res_penalty_terms',...
%     'valid_all_penalty_terms','valid_all_sampled_params'); 
% valid_high_res_penalty_terms = result_file.valid_high_res_penalty_terms; 
% valid_all_penalty_terms = result_file.valid_all_penalty_terms; 
% valid_all_sampled_params = result_file.valid_all_sampled_params; 
% 
% rank_valid_all_penalty_terms = tiedrank(valid_high_res_penalty_terms);
% rank_param_mat = tiedrank(valid_all_sampled_params); 
% augmented_PRCC = partialcorri(rank_valid_all_penalty_terms,rank_param_mat); 
% 
%     % Create problemObject
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject_ctrl = setProblemObject_v3(group_description,[],param_info_path); 
% estimated_param_name_list = {problemObject_ctrl.Estimated.Name}; 
% 
%     % Filter out parameters that are sensitive to at least 1 penalty terms 
% param_relax_idx_list = [];
% param_relax_name_list = {}; 
% for param_idx = 1:size(augmented_PRCC,2)
%     param_corr = augmented_PRCC(:,param_idx);
%     if any(abs(param_corr) > 0.15)
%         param_relax_idx_list = [param_relax_idx_list param_idx];
%         param_relax_name_list{end + 1} = estimated_param_name_list{param_idx}; 
%     end
% end
% 
%     % Relax their bounds 
%         % Get original bounds from problemObject_ctrl
% param_bounds = [problemObject_ctrl.Estimated.Bounds]; 
% param_bounds = reshape(param_bounds,[2,length(estimated_param_name_list)]);
% param_bounds = param_bounds'; 
% relaxed_param_bounds = param_bounds; 
% for relax_idx = 1:length(param_relax_idx_list)
% 
%     param_relax_idx = param_relax_idx_list(relax_idx); 
%     relaxed_param_bounds(relax_idx,:) = [param_bounds(relax_idx,1) ./ 1e+03,param_bounds(relax_idx,2) .* 1e+03]; % Expand 3 orders-of-magnitude up and down  
% 
% end
% 
%     % Redo LHS sampling
% num_sample = 1e+06;
% save_path = 'test_save_files/20240822_lhs_sampled_params_sensitivity_exp.mat';
% lhs_sample(problemObject_ctrl,relaxed_param_bounds,num_sample,save_path); 
% 
% % Try optimization-based bound relaxation 
% % lhs_sample(problemObject_ctrl,param_bounds,3000,'test_save_files/20240822_lhs_sampled_params_optimization.mat'); 

%% Parameter + initial condition sampling with expanded range for parameters 

% % Modify the bounds in parameters_v2.xlsx to 0.1x for LB and 10x for UB 
% param_info_table = readtable('param_info/parameters_v2.xlsx'); 
% for row_idx = 1:height(param_info_table)
% 
%     param_info_table.LB(row_idx) = 0.01 .* param_info_table.LB(row_idx); 
%     param_info_table.UB(row_idx) = 100 .* param_info_table.UB(row_idx); 
% 
% end
% 
% writetable(param_info_table,'param_info/parameters_v2_expanded.xlsx'); 
% 
% % Sample using the updated parameter ranges 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2_expanded.xlsx';
% problemObject = setProblemObject_v3(group_description,[],param_info_path); 
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% modelObj = problemObject.Model; 
% experimental_grouped_data = problemObject.Data; 
% % simFunction = create_simFun_from_problemObject(problemObject);
% % dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
% num_kinetic_params = length(problemObject.Estimated); 
% 
% kinetic_param_names = {problemObject.Estimated.Name};
% kinetic_param_bounds = [problemObject.Estimated.Bounds]; 
% kinetic_param_bounds = reshape(kinetic_param_bounds,[2,num_kinetic_params]); 
% kinetic_param_bounds = kinetic_param_bounds'; 
% 
% % % Fix mRNA degradation parameters (using 95% CI as LB/UB) 
% % mRNA_deg_param_info_file = load('test_save_files/20240530_mRNA_degradation_param_bounds.mat'); 
% % for estimated_param_idx = 1:length(mRNA_deg_param_info_file.estimated_params_labels)
% %     estimated_param_label = mRNA_deg_param_info_file.estimated_params_labels{estimated_param_idx};
% %     estimated_param_bounds = mRNA_deg_param_info_file.param_bounds_CI(estimated_param_idx,:);
% % 
% %     kinetic_param_idx = strcmp(estimated_param_label,kinetic_param_names); 
% %     kinetic_param_bounds(kinetic_param_idx,:) = estimated_param_bounds; 
% % 
% % end
% log_param_bounds = log10(kinetic_param_bounds); 
% 
% % Generate sampled kinetic parameter values 
% num_sample = 1000000; % 1e+06 
% lhsSamples = lhsdesign(num_sample,num_kinetic_params); 
% sampled_params = nan(num_sample,num_kinetic_params); 
% for i = 1:num_kinetic_params
%     lowerBound = log_param_bounds(i, 1);
%     upperBound = log_param_bounds(i, 2);
%     sampled_params(:, i) = lowerBound + (upperBound - lowerBound) * lhsSamples(:, i);
% end
% 
% nominal_sampled_params = 10 .^ sampled_params;
% 
%     % Define parameters to be constrained/filtered 
% kinetic_param_name_oi_list = {'TXTL_RNAdeg_kc_sfGFP','TXTL_RNAdeg_kc_kanR','TXTL_RNAdeg_kc',...
%     'TXTL_RNAdeg_R_sfGFP','TXTL_RNAdeg_R','TXTL_PT7_RNAPbound_R','TXTL_PT773_RNAPbound_R',...
%     'TXTL_PJ23119_RNAPbound_R','TXTL_PJ23105_RNAPbound_R'}; 
%     % Find their indices in params list 
% kinetic_param_name_oi_idx_list = nan(length(kinetic_param_name_oi_list),1); 
% for kinetic_param_oi_idx = 1:length(kinetic_param_name_oi_idx_list)
%     kinetic_param_name_oi = kinetic_param_name_oi_list{kinetic_param_oi_idx}; 
%     kinetic_param_name_oi_idx_list(kinetic_param_oi_idx) = find(strcmp(kinetic_param_name_oi,kinetic_param_names)); 
% end
% 
%     % Encode a matrix to represent relationships among parameters
%         % Each row represents a criteria and each column represent
%         % parameters involved 
% criteria_matrix = [1  0  -1  0  0  0  0  0  0; 
%                    0  1  -1  0  0  0  0  0  0;
%                    0  0   0 -1  1  0  0  0  0;
%                    0  0   0  0  0  1 -1  0  0;
%                    0  0   0  0  0  0  1  0  -1;
%                    0  0   0  0  0  0  0  1  -1;]; % Assumes < 0  
% satisfied_idx_list = []; % Initialize 
% for criterion_idx = 1:size(criteria_matrix,1)
%     criteria_vec = criteria_matrix(criterion_idx,:); 
%     smaller_param_idxa = find(criteria_vec == 1); 
%     larger_param_idxa = find(criteria_vec == -1); 
%     row_satisfied_idx_list = nominal_sampled_params(:,kinetic_param_name_oi_idx_list(smaller_param_idxa)) < nominal_sampled_params(:,kinetic_param_name_oi_idx_list(larger_param_idxa));
%     if isempty(satisfied_idx_list)
%         satisfied_idx_list = row_satisfied_idx_list; 
%     else
%         satisfied_idx_list = satisfied_idx_list & row_satisfied_idx_list; 
%     end
% 
% end
% final_satisfied_idx_list = find(satisfied_idx_list);
% 
% satisfied_sampled_params = nominal_sampled_params(final_satisfied_idx_list,:); 
% 
% save('test_save_files/20240819_lhs_sampled_params.mat','satisfied_sampled_params','problemObject','dosing_information'); 
% for param_idx = 1:1 %size(satisfied_sampled_params,1)
% 
%     sampled_params_selected = satisfied_sampled_params(param_idx,:); 
% 
%     [scaled_penalty_vec,penalty_term_labels,penalty_term_length,unscaled_penalty_vec] = wrapper_calculate_obj_higher_resolution(sampled_params_selected,...
%         problemObject,dosing_information,[],[],mode,num_conc,num_promotor);
% end