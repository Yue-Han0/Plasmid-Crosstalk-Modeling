% clear
% clc

currentpath = pwd; 
addpath(genpath(currentpath))
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))
promotor_name_list = {'T7_strong','T7_weak','sigma70_strong','sigma70_weak'};
conc_vec = [0.5,1,2.5,5,10,15,30]; 
mode = 'PE';
num_conc = 7;
num_promotor = 4; 
num_comb = 4; 
all_obj_string_file = load('qualitative_trend_obj_string.mat'); 
all_obj_string = all_obj_string_file.all_obj_string; 
all_obj_string = all_obj_string'; 
tStart = 0; 
tEnd = 21600; % Modify this later to draw info from data  

%% Compile results for Conflict in Achieving Large Positive Crosstalk in Weak Promotor
%     % Evidence # 1: single objective optimization 
% directory_path = 'param_est_run_save'; 
% directory_files = dir(directory_path); 
% file_names = {directory_files.name};
% % for obj_idx = 16:22 % single objective indices 
% %     resultFileName_idx = contains(file_names,sprintf('objString_%d',obj_idx)); 
% % 
% % end
% sample_obj_string_16_result_file = load('param_est_run_save/20240516_param_est_run1535_objString_16.mat');
% 
%     % Evidence # 2: Goal Attainment Optimization 
% % Define a list of goals and weights 
%         % penalty_vec = (1) baseline_data_dev (2) crosstalk_ratio_dev (3)
%     % T7strong toxin (4) T7weak toxin (5) sigma70strong toxin (6)
%     % sigma70weak no toxin (7) positive crosstalk (8) negative
%     % crosstalk (9) promotor strength (apply to positive crosstalk
%     % only) (10) other penalties (11) Promotor sensitivity 
% % Hardcode this for now to quickly try this out 
% % goal = [10;2;-2;-2;-2;-2;-0.2 .* ones(12,1);-0.2 .* ones(12,1);-2;0;-10];
% goal_attainment_sampled_result_file = load('param_est_run_save/20240531_param_est_run0918_goal_attainment_run1_iter_2.mat'); 
% 
%     % Evidence # 3: Large-scale Sampling 
% large_scale_sampling_result_file = load('test_save_files/20240619_large_scale_parameter_sampling_summary.mat'); 
% all_penalty_terms = large_scale_sampling_result_file.all_penalty_terms; 


%% Compile sampling results; Case study on nan values; etc., 
    
%     % Check for missing files 
% for iter = 1:1000
%     result_file_name = sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d.mat',iter); 
%     if ~exist(result_file_name,'file')
%         fprintf('\n Missing run iteration #%d',iter)
%     end
% end
% 
% for iter = 11:1000
%     mislabeled_result_file_name = sprintf('param_est_run_save/20240717_param_and_init_cond_sampling_constrained_run%d.mat',iter); 
%     if exist(mislabeled_result_file_name,'file')
%         actual_result_file_name = sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d.mat',iter); 
%         if ~exist(actual_result_file_name,'file')
%             movefile(mislabeled_result_file_name,actual_result_file_name)
%         end
%     end
% end
%     % Test code to compile results 
%     % Load a sample result file 
% sample_result_file_name = 'param_est_run_save/20240717_param_and_init_cond_sampling_constrained_run1.mat'; 
% sample_result_file = load(sample_result_file_name); 
% num_params = size(sample_result_file.sampled_params_selected,2); 
% num_penalty_terms = size(sample_result_file.all_penalty_terms,2); 
% num_high_res_penalty_terms = size(sample_result_file.high_res_all_penalty_terms,2); 
% penalty_labels = sample_result_file.penalty_labels; 
% high_res_penalty_labels = sample_result_file.high_res_penalty_term_labels; 
% high_res_penalty_length = sample_result_file.high_res_penalty_term_length; 
% 
% num_params_sampled = 10000; 
% num_init_cond_sampled = 100; 
% 
%     % Preassign 
% all_sampled_params = nan(num_params_sampled * num_init_cond_sampled,num_params); 
% all_calculated_penalty_terms = nan(num_params_sampled * num_init_cond_sampled,num_penalty_terms); 
% all_high_res_penalty_terms = nan(num_params_sampled * num_init_cond_sampled,num_high_res_penalty_terms); 
% all_high_res_penalty_terms_unscaled = nan(num_params_sampled * num_init_cond_sampled,num_high_res_penalty_terms); 
% 
% for iter = 1:1000
%         % Load result file 
%     try 
%         result_file_name = sprintf('param_est_run_save/20240717_param_and_init_cond_sampling_constrained_run%d.mat',iter); 
%         result_file = load(result_file_name); 
%         sampled_params_selected = result_file.sampled_params_selected;
%         all_penalty_terms = result_file.all_penalty_terms; 
%         high_res_penalty_terms = result_file.high_res_all_penalty_terms; 
%         high_res_penalty_terms_unscaled = result_file.high_res_all_penalty_terms_unscaled;
%     catch
%         fprintf('\n Missing Run # %d',iter); 
%     end
% 
%         % Assign 
%     for local_param_iter = 1:size(sampled_params_selected,1)
% 
%         start_idx = (iter - 1) * 10 * num_init_cond_sampled + (local_param_iter - 1) * num_init_cond_sampled + 1; 
%         stop_idx = (iter - 1) * 10 * num_init_cond_sampled + local_param_iter * num_init_cond_sampled; 
% 
%         sampled_params_single = sampled_params_selected(local_param_iter,:); 
%         all_sampled_params(start_idx:stop_idx,:) = repmat(sampled_params_single,num_init_cond_sampled,1);
%         all_calculated_penalty_terms(start_idx:stop_idx,:) = all_penalty_terms((local_param_iter - 1) * num_init_cond_sampled + 1:local_param_iter * num_init_cond_sampled,:);
%         all_high_res_penalty_terms(start_idx:stop_idx,:) = high_res_penalty_terms((local_param_iter - 1) * num_init_cond_sampled + 1:local_param_iter * num_init_cond_sampled,:);
%         all_high_res_penalty_terms_unscaled(start_idx:stop_idx,:) = high_res_penalty_terms_unscaled((local_param_iter - 1) * num_init_cond_sampled + 1:local_param_iter * num_init_cond_sampled,:);
%     end
% 
% end
% 
% save('test_save_files/20240717_sampling_constrained_result_summary.mat','all_calculated_penalty_terms','all_sampled_params','penalty_labels',...
%     'all_high_res_penalty_terms','all_high_res_penalty_terms_unscaled','high_res_penalty_labels','high_res_penalty_length');

%     % Case study for cause of nan values in penalty terms - numerical
%     % instability? 
% sampled_params_file = load('test_save_files/20240613_param_init_cond_sampling_constrained.mat'); 
% rerun_flag_result_file = load('test_save_files/20240715_samping_missing_file_check.mat'); 
% nan_containing_idx_list = find(rerun_flag_result_file.run_flag == 2); 
% 
% case_study_1_file = load(sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d.mat',5));
% case_study_2_file = load(sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d.mat',7));
%     % For these two case studies,  nan values found for certain sampled
%     % parameters, find those parameters 
% nan_row_idx_list = find(any(isnan(case_study_1_file.all_penalty_terms),2)); 
% corresponding_param_idx = unique(ceil(nan_row_idx_list ./ 100)); 
% sampled_params = case_study_1_file.sampled_params_selected(corresponding_param_idx,:); 
%     % Create simFunction 
% model = sampled_params_file.modelObj; 
% sampled_param_names = sampled_params_file.kinetic_param_names; 
% species_name_list = {model.Species.Name}; 
% remove_idx = strcmp('protein sfGFP*',species_name_list); 
% additional_track_species_list = species_name_list(~remove_idx); 
% observables = [{'protein sfGFP*'},additional_track_species_list];
% simFunction = createSimFunction(model,sampled_param_names,observables,sampled_params_file.target_name_list); 
%     % Extract dosing information and run simulation 
% dosing_information = case_study_1_file.all_dosing_information{1}; % choose first one
% [simulated_time,simulated_data] = simFunction(sampled_params,tEnd,dosing_information,tStart:tEnd);


% %     % Update stored model in saved sampled paraemters, and rerun sampling using the updated model (without duplicate reactions) 
% sampled_params_file = load('20240613_param_init_cond_sampling_constrained.mat'); 
%     % Get necessary variables from sampled params file 
% nominal_sampled_params = sampled_params_file.nominal_sampled_params; 
% satisfied_idx_list = sampled_params_file.satisfied_idx_list; 
% all_dosing_information = sampled_params_file.all_dosing_information; 
% experimental_grouped_data = sampled_params_file.experimental_grouped_data; 
% kinetic_param_names = sampled_params_file.kinetic_param_names; 
% target_name_list = sampled_params_file.target_name_list; 
%     % Get the new model 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject = setProblemObject_v3(group_description,[],param_info_path);
% modelObj = problemObject.Model; 
%     % Save sampled params and initial conditions to new file 
% save('test_save_files/20240717_param_init_cond_sampling_constrained.mat'); 


%% Check updated sampling files 
% run_flag = nan(1000,1); 
% for iter = 1:1000
% 
%     result_file_name = sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d.mat',iter); 
%     try
%         result_file = load(result_file_name); 
%         if ~any(any(isnan(result_file.all_penalty_terms)))
%             run_flag(iter) = 1; % all good 
%         else
%             run_flag(iter) = 2; % calculated penlaty value containing nan 
%         end
%     catch
%         run_flag(iter) = 0; % result file missing 
%     end
% 
% end

% save('test_save_files/20240715_samping_missing_file_check.mat','run_flag'); 

% load('test_save_files/20240715_samping_missing_file_check.mat','run_flag');

% If result file missing, resubmit job; if all penalty value contains nan,
% investigate more 
% command = fopen('jobSubmission/sample_params_and_init_cond_account_round2.txt','w');
% for run_idx = 1:length(run_flag)
%     % rerun_flag = true; 
%     % result_file_name = sprintf('param_est_run_save/20240613_param_and_init_cond_sampling_constrained_run%d.mat',run_idx);
%     % if exist(result_file_name,'file')
%     %     rerun_flag = false; 
%     % end
%     if isequal(run_flag(run_idx),0)
%         fprintf(command,'\n sbatch -A gts-mstyczynski6 -q inferno -N1 --ntasks-per-node=1 --mem-per-cpu=16G -t 72:00:00 jobSubmission/sample_params_and_init_cond_run_%d.txt',run_idx);
%     end
% end
% for run_idx = 1:length(run_flag)
%     if isequal(run_flag(run_idx),2)
%         % Load file and see if all calculated penalty values are nan?
%         result_file_name = sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d.mat',run_idx); 
%         result_file = load(result_file_name);
% 
%     end
% end

%% Large Crosstalk Ratio Incompatibility Check 

% Experiment #1: Compare fitted parameters to lower plasmid
% concentration and the large-scale sampling / optimization 

%     % Load result files 
% large_scale_sampling_result_summary_file = load('test_save_files/20240619_large_scale_parameter_sampling_summary.mat'); 
% large_scale_fitting_result_summary_file = load('param_est_run_save/20231106_large_scale_run.mat'); 
% RFM_integration_ctrl_result_file = load('param_est_run_save/20230709_param_est_run1945_RFM_integration_ctrl.mat'); 
% fitting_result_summary = large_scale_fitting_result_summary_file.updated_optimization_run_table; 
% all_sampled_params = large_scale_sampling_result_summary_file.all_sampled_params; 
% all_penalty_terms = large_scale_sampling_result_summary_file.all_penalty_terms; 

%     % Check RFM control results and only select sets that reproduce large
%     % positive crosstalk (Ratio > 2.5) - Best 1-7 reproduces large positive
%     % crosstalk ratio 
% [~,sorted_SSE_idx_ctrl] = sort(RFM_integration_ctrl_result_file.metric_summary,'ascend'); 
% species_oi = 'protein sfGFP*'; 
% conc_vec_oi = conc_vec(1:3); 
% % Plot out the experimental data 
% problemObject_ctrl = RFM_integration_ctrl_result_file.problemObject;
% experimental_groupedData = problemObject_ctrl.Data; 
% 
% for idx = 1:48
%     % Select fitResult of interest and get simData
%     fitResult_oi_ctrl = RFM_integration_ctrl_result_file.all_fitResults{sorted_SSE_idx_ctrl(idx)}; 
%     simData_ctrl = fitted(fitResult_oi_ctrl);
% 
%     % These should be organized as no empty - conc_vec 1-3, empty -
%     % conc_vec 1-3 ? Let's try this first and plot out the time-course for
%     % GFP 
%     figure;
%     for conc_idx = 1:length(conc_vec_oi)
% 
%         subplot(2,2,conc_idx)
% 
%         simData_no_empty_ctrl = simData_ctrl(conc_idx); 
%         simData_empty_ctrl = simData_ctrl(length(conc_vec_oi) + conc_idx); 
%         species_oi_idx_ctrl = find(strcmp(simData_no_empty_ctrl.DataNames,species_oi));
% 
%         % Plot out experimental data 
%         group_num_no_empty = 168 + conc_idx; 
%         group_num_empty = 168 + 7 + conc_idx; 
%         no_empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_no_empty); 
%         no_empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_no_empty); 
%         empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_empty); 
%         empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_empty); 
%         plot(no_empty_timeVec,no_empty_GFP_conc,'--','Color','k','LineWidth',1.5)
%         hold on
%         plot(empty_timeVec,empty_GFP_conc,'-','Color','k','LineWidth',1.5)
% 
% 
%         plot(simData_no_empty_ctrl.Time,simData_no_empty_ctrl.Data(:,species_oi_idx_ctrl),'--','LineWidth',1.5,'Color','g')
%         plot(simData_empty_ctrl.Time,simData_empty_ctrl.Data(:,species_oi_idx_ctrl),'-','LineWidth',1.5,'Color','g')
% 
%         title(sprintf('%.1f nM Reporter Plasmid',conc_vec_oi(conc_idx)))
% 
%     end
%     legend('Exp - no empty','Exp - empty','Ctrl - no empty','Ctrl - empty')
%     sgtitle(sprintf('Best Fitting #%d',idx))
% end

%     % Plot out parameter distribution for all
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% 
%     % List of parameter names for sampling
% problemObject_sampling = setProblemObject_v2(group_description,param_info_path,[]); 
% estimated_param_name_list_sampling = {problemObject_sampling.Estimated.Name};
% 
%     % List of parameter names for fitting 
% % sample_fitting_result_file = load('param_est_run_save/20231102_param_est_run1531_1.mat'); 
% problemObject_fitting = sample_fitting_result_file.problemObject; 
% estimated_param_name_list_fitting = {problemObject_fitting.Estimated.Name}; 
% 
%     % Compile estimated parameter names for weak promotor low
%     % concentrations 
% sigma70_weak_param_name_list = {RFM_integration_ctrl_result_file.problemObject.Estimated.Name};
% sigma70_weak_estimated_params = nan(length(RFM_integration_ctrl_result_file.metric_summary),length(sigma70_weak_param_name_list)); 
% for iter = 1:length(RFM_integration_ctrl_result_file.metric_summary)
% 
%     fitResult_oi = RFM_integration_ctrl_result_file.all_fitResults{iter}; 
%     sigma70_weak_estimated_params_single = [fitResult_oi.ParameterEstimates.Estimate]; 
%     sigma70_weak_estimated_params(iter,:) = sigma70_weak_estimated_params_single; 
% 
% end
%     % get sorted metric summary 
% [~,sorted_SSE_idx_ctrl] = sort(RFM_integration_ctrl_result_file.metric_summary,'ascend'); 
% 
%     % 
% figure; 
% for param_idx = 1:length(estimated_param_name_list_fitting)
% 
%     subplot(ceil(sqrt(length(estimated_param_name_list_fitting))),ceil(sqrt(length(estimated_param_name_list_fitting))),param_idx)
% 
%     param_name_oi = estimated_param_name_list_fitting{param_idx}; 
%     param_name_idx_sampling = find(strcmp(param_name_oi,estimated_param_name_list_sampling)); 
%     param_name_idx_sigma70_weak = find(strcmp(param_name_oi,sigma70_weak_param_name_list)); 
%     sampled_params = all_sampled_params(:,param_name_idx_sampling); 
%     estimated_params = nan(height(fitting_result_summary),1); 
%     for row_idx = 1:height(fitting_result_summary)
%         estimated_params_all = fitting_result_summary.ParameterEstimates(row_idx); 
%         estimated_params_all = estimated_params_all{1,1}; 
%         estimated_params(row_idx) = estimated_params_all(param_idx); 
%     end
% 
%     histogram(log10(estimated_params))
%     hold on 
%     % if ~isempty(param_name_idx_sampling)
%     %     histogram(log10(sampled_params))
%     % end
% 
%     % plot estimated value in weak promotor low conc fitting 
%     if ~isempty(param_name_idx_sigma70_weak)
%         for j = 1:7
%             sigma70_weak_idx_oi = sorted_SSE_idx_ctrl(j); 
%             sigma70_weak_estimated_params_single = sigma70_weak_estimated_params(sigma70_weak_idx_oi,:); 
%             xline(log10(sigma70_weak_estimated_params_single(param_name_idx_sigma70_weak)),'LineWidth',1.5)
%         end
%     end
%     title(strrep(param_name_oi,'_',' '))
% 
% end

%% Check RFM integration results round #2 
% RFM_integration_ctrl_result_file = load('param_est_run_save/20230709_param_est_run1945_RFM_integration_ctrl.mat'); 
% RFM_integration_test_result_file = load('param_est_run_save/20230709_param_est_run1528_RFM_integration.mat'); 

%     % Compile fit results 
% wrapper_analyze_fit_result_v2(RFM_integration_ctrl_result_file.all_fitResults,RFM_integration_ctrl_result_file.problemObject,'param_est_run_save/20230709_param_est_run1945_RFM_integration_ctrl.mat','SSE',false,false) 
% wrapper_analyze_fit_result_v2(RFM_integration_test_result_file.all_fitResults,RFM_integration_test_result_file.problemObject,'param_est_run_save/20230709_param_est_run1528_RFM_integration.mat','SSE',false,false) 

%     % Compare fitted SSE distribution & lowest SSE 
% figure;
% histogram(log10(RFM_integration_ctrl_result_file.metric_summary))
% hold on 
% histogram(log10(RFM_integration_test_result_file.metric_summary))
% legend('Control','RFM integration')
% 
% fprintf('\n MinSSE Control: %.2f',min(RFM_integration_ctrl_result_file.metric_summary))
% fprintf('\n MinSSE Test: %.2f',min(RFM_integration_test_result_file.metric_summary))
% 
% 
% % Plot out the experimental data 
% problemObject_ctrl = RFM_integration_ctrl_result_file.problemObject;
% % problemObject_test = RFM_integration_test_result_file.problemObject; 
% experimental_groupedData = problemObject_ctrl.Data; 
%     % For each result file, select top 5 best fitted model and check crosstalk
%     % ratio 
% [~,sorted_SSE_idx_ctrl] = sort(RFM_integration_ctrl_result_file.metric_summary,'ascend'); 
% [~,sorted_SSE_idx_test] = sort(RFM_integration_test_result_file.metric_summary,'ascend'); 
% 
% species_oi = 'protein sfGFP*'; 
% conc_vec_oi = conc_vec(1:3); 
% 
% for idx = 1:5
%     % Select fitResult of interest and get simData
%     fitResult_oi_ctrl = RFM_integration_ctrl_result_file.all_fitResults{sorted_SSE_idx_ctrl(idx)}; 
%     fitResult_oi_test = RFM_integration_test_result_file.all_fitResults{sorted_SSE_idx_test(idx)}; 
%     simData_ctrl = fitted(fitResult_oi_ctrl);
%     simData_test = fitted(fitResult_oi_test); 
% 
%     % These should be organized as no empty - conc_vec 1-3, empty -
%     % conc_vec 1-3 ? Let's try this first and plot out the time-course for
%     % GFP 
%     figure;
%     for conc_idx = 1:length(conc_vec_oi)
% 
%         subplot(2,2,conc_idx)
% 
%         simData_no_empty_ctrl = simData_ctrl(conc_idx); 
%         simData_no_empty_test = simData_test(conc_idx);
%         simData_empty_ctrl = simData_ctrl(length(conc_vec_oi) + conc_idx); 
%         simData_empty_test = simData_test(length(conc_vec_oi) + conc_idx);
% 
%         species_oi_idx_ctrl = find(strcmp(simData_no_empty_ctrl.DataNames,species_oi));
%         species_oi_idx_test = find(strcmp(simData_no_empty_test.DataNames,species_oi));
% 
%         % Plot out experimental data 
%         group_num_no_empty = 168 + conc_idx; 
%         group_num_empty = 168 + 7 + conc_idx; 
%         no_empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_no_empty); 
%         no_empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_no_empty); 
%         empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_empty); 
%         empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_empty); 
%         plot(no_empty_timeVec,no_empty_GFP_conc,'--','Color','k','LineWidth',1.5)
%         hold on
%         plot(empty_timeVec,empty_GFP_conc,'-','Color','k','LineWidth',1.5)
% 
% 
%         plot(simData_no_empty_ctrl.Time,simData_no_empty_ctrl.Data(:,species_oi_idx_ctrl),'--','LineWidth',1.5,'Color','g')
%         hold on 
%         plot(simData_no_empty_test.Time,simData_no_empty_test.Data(:,species_oi_idx_test),'--','LineWidth',1.5,'Color','r')
%         plot(simData_empty_ctrl.Time,simData_empty_ctrl.Data(:,species_oi_idx_ctrl),'-','LineWidth',1.5,'Color','g')
%         plot(simData_empty_test.Time,simData_empty_test.Data(:,species_oi_idx_test),'-','LineWidth',1.5,'Color','r')
% 
%         title(sprintf('%.1f nM Reporter Plasmid',conc_vec_oi(conc_idx)))
% 
%     end
%     legend('Exp - no empty','Exp - empty','Ctrl - no empty','Test - no empty','Ctrl - empty','Test - empty')
%     sgtitle(sprintf('Best Fitting #%d',idx))
% end

% % Plot out more species in cases where large crosstalk ratio is captured 
% [~,sorted_SSE_idx_ctrl] = sort(RFM_integration_ctrl_result_file.metric_summary,'ascend'); 
% [~,sorted_SSE_idx_test] = sort(RFM_integration_test_result_file.metric_summary,'ascend'); 
% 
%     % Resource of interest: RNAP, Ribo, RNase, AGTP, CUTP, AA, RNA
%     % utrGFP-sfGFP (not bound to RNase), RNA sfGFP bound to RNase,
%     % RNA utrkanR-kanR (not bound to RNase), RNA kanR bound to
%     % RNase, protein kanR, protein sfGFP*
% resource_plot_names = {'RNAP','Ribo','RNase','AGTP','CUTP','AA','protein kanR','protein sfGFP*','toxin','RNA sfGFP not bound to RNase',...
%     'RNA sfGFP bound to RNase','RNA kanR not bound to RNase','RNA kanR bound to RNase'};
% all_species_name_ctrl = {RFM_integration_ctrl_result_file.problemObject.Model.Species.Name}; 
% all_species_name_test = {RFM_integration_test_result_file.problemObject.Model.Species.Name}; 
% 
%     % Compile species list for sum 
% RNA_sfGFP_not_bound_idx_list_ctrl = contains(all_species_name_ctrl,'RNA utrGFP--sfGFP') & ~contains(all_species_name_ctrl,'RNase'); 
% RNA_sfGFP_bound_idx_list_ctrl = contains(all_species_name_ctrl,'RNA utrGFP--sfGFP') & contains(all_species_name_ctrl,'RNase'); 
% RNA_kanR_not_bound_idx_list_ctrl = contains(all_species_name_ctrl,'RNA utrkanR--kanR') & ~contains(all_species_name_ctrl,'RNase'); 
% RNA_kanR_bound_idx_list_ctrl = contains(all_species_name_ctrl,'RNA utrkanR--kanR') & contains(all_species_name_ctrl,'RNase'); 
% 
% RNA_sfGFP_not_bound_idx_list_test = contains(all_species_name_test,'RNA utrGFP--sfGFP') & ~contains(all_species_name_test,'RNase'); 
% RNA_sfGFP_bound_idx_list_test = contains(all_species_name_test,'RNA utrGFP--sfGFP') & contains(all_species_name_test,'RNase'); 
% RNA_kanR_not_bound_idx_list_test = contains(all_species_name_test,'RNA utrkanR--kanR') & ~contains(all_species_name_test,'RNase'); 
% RNA_kanR_bound_idx_list_test = contains(all_species_name_test,'RNA utrkanR--kanR') & contains(all_species_name_test,'RNase'); 
% 
% RNA_sfGFP_not_bound_species_ctrl = all_species_name_ctrl(RNA_sfGFP_not_bound_idx_list_ctrl); 
% RNA_sfGFP_bound_species_ctrl = all_species_name_ctrl(RNA_sfGFP_bound_idx_list_ctrl); 
% RNA_kanR_not_bound_species_ctrl = all_species_name_ctrl(RNA_kanR_not_bound_idx_list_ctrl); 
% RNA_kanR_bound_species_ctrl = all_species_name_ctrl(RNA_kanR_bound_idx_list_ctrl); 
% 
% RNA_sfGFP_not_bound_species_test = all_species_name_test(RNA_sfGFP_not_bound_idx_list_test); 
% RNA_sfGFP_bound_species_test = all_species_name_test(RNA_sfGFP_bound_idx_list_test); 
% RNA_kanR_not_bound_species_test = all_species_name_test(RNA_kanR_not_bound_idx_list_test); 
% RNA_kanR_bound_species_test = all_species_name_test(RNA_kanR_bound_idx_list_test); 
% 
% resource_model_names_ctrl = {'RNAP','Ribo','RNase','AGTP','CUTP','AA','protein kanR','protein sfGFP*','toxin',{RNA_sfGFP_not_bound_species_ctrl},...
%     {RNA_sfGFP_bound_species_ctrl},{RNA_kanR_not_bound_species_ctrl},{RNA_kanR_bound_species_ctrl}}; 
% resource_model_names_test = {'RNAP','Ribo','RNase','AGTP','CUTP','AA','protein kanR','protein sfGFP*','toxin',{RNA_sfGFP_not_bound_species_test},...
%     {RNA_sfGFP_bound_species_test},{RNA_kanR_not_bound_species_test},{RNA_kanR_bound_species_test}}; 
% 
% for idx = 1:1
% 
%     % Select fitResult of interest and get simData
%         % RFM Control 
%     fitResult_oi_ctrl = RFM_integration_ctrl_result_file.all_fitResults{sorted_SSE_idx_ctrl(idx)}; 
%     simData_ctrl = fitted(fitResult_oi_ctrl);
% 
%     simData_ctrl_no_empty_conc_vec1 = simData_ctrl(1); 
%     simData_ctrl_empty_conc_vec1 = simData_ctrl(4); 
% 
%         % RFM test 
%     fitResult_oi_test = RFM_integration_test_result_file.all_fitResults{sorted_SSE_idx_test(idx)}; 
%     simData_test = fitted(fitResult_oi_test);
% 
%     simData_test_no_empty_conc_vec1 = simData_test(1); 
%     simData_test_empty_conc_vec1 = simData_test(4); 
% 
%     figure;
%         % Plot out all species 
%     for species_idx = 1:length(simData_test_no_empty_conc_vec1.DataNames)
% 
%         subplot(ceil(sqrt(length(simData_test_no_empty_conc_vec1.DataNames))),ceil(sqrt(length(simData_test_no_empty_conc_vec1.DataNames))),species_idx)
% 
%         species_name = simData_test_no_empty_conc_vec1.DataNames{species_idx}; 
%         plot(simData_test_no_empty_conc_vec1.Time,simData_test_no_empty_conc_vec1.Data(:,species_idx),'LineWidth',1.5,'Color','g')
%         hold on 
%         plot(simData_test_empty_conc_vec1.Time,simData_test_empty_conc_vec1.Data(:,species_idx),'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         title(strrep(species_name,'_',' '))
% 
%     end
% 
%     %     % Plot all resources of interest 
%     % for resource_idx = 1:length(resource_plot_names)
%     % 
%     %     subplot(ceil(sqrt(length(resource_plot_names))),ceil(sqrt(length(resource_plot_names))),resource_idx)
%     % 
%     %     resource_plot_name = resource_plot_names{resource_idx}; 
%     %     resource_model_name_ctrl = resource_model_names_ctrl{resource_idx}; 
%     %     resource_model_name_test = resource_model_names_test{resource_idx}; 
%     %     % if ischar(resource_model_name_ctrl)
%     %     %     data_idx_ctrl = strcmp(simData_ctrl_no_empty_conc_vec1.DataNames,resource_model_name_ctrl); 
%     %     %     plot(simData_ctrl_no_empty_conc_vec1.Time,simData_ctrl_no_empty_conc_vec1.Data(:,data_idx_ctrl),'LineWidth',1.5,'Color','g')
%     %     %     hold on 
%     %     %     plot(simData_ctrl_empty_conc_vec1.Time,simData_ctrl_empty_conc_vec1.Data(:,data_idx_ctrl),'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%     %     % else
%     %     %     resource_oi_model_name_list_ctrl = resource_model_name_ctrl{1,1}; 
%     %     %     data_idx_list_ctrl = nan(length(resource_oi_model_name_list_ctrl),1);
%     %     %     for resource_model_idx = 1:length(resource_oi_model_name_list_ctrl)
%     %     %         resource_oi_model_name = resource_oi_model_name_list_ctrl{resource_model_idx}; 
%     %     %         data_idx_list_ctrl(resource_model_idx) = find(strcmp(resource_oi_model_name,simData_ctrl_no_empty_conc_vec1.DataNames)); 
%     %     %     end
%     %     %     plot(simData_ctrl_no_empty_conc_vec1.Time,sum(simData_ctrl_no_empty_conc_vec1.Data(:,data_idx_list_ctrl),2),'LineWidth',1.5,'Color','g')
%     %     %     hold on 
%     %     %     plot(simData_ctrl_empty_conc_vec1.Time,sum(simData_ctrl_empty_conc_vec1.Data(:,data_idx_list_ctrl),2),'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%     %     % 
%     %     % end
%     %     if ischar(resource_model_name_test)
%     %         data_idx_test = strcmp(simData_test_no_empty_conc_vec1.DataNames,resource_model_name_test); 
%     %         plot(simData_test_no_empty_conc_vec1.Time,simData_test_no_empty_conc_vec1.Data(:,data_idx_test),'LineWidth',1.5,'Color','g','LineStyle','--')
%     %         plot(simData_test_empty_conc_vec1.Time,simData_test_empty_conc_vec1.Data(:,data_idx_test),'LineWidth',1.5,'Color',[0.5,0.5,0.5],'LineStyle','--')
%     %     else
%     %         resource_oi_model_name_list_test = resource_model_name_test{1,1}; 
%     %         data_idx_list_test = nan(length(resource_oi_model_name_list_test),1);
%     %         for resource_model_idx = 1:length(resource_oi_model_name_list_test)
%     %             resource_oi_model_name = resource_oi_model_name_list_test{resource_model_idx}; 
%     %             data_idx_list_test(resource_model_idx) = find(strcmp(resource_oi_model_name,simData_test_no_empty_conc_vec1.DataNames)); 
%     % 
%     %         end
%     %         plot(simData_test_no_empty_conc_vec1.Time,sum(simData_test_no_empty_conc_vec1.Data(:,data_idx_list_test),2),'LineWidth',1.5,'Color','g','LineStyle','--')
%     %         plot(simData_test_empty_conc_vec1.Time,sum(simData_test_empty_conc_vec1.Data(:,data_idx_list_test),2),'LineWidth',1.5,'Color',[0.5,0.5,0.5],'LineStyle','--')
%     %     end
%     % 
%     %     title(strrep(resource_plot_name,'_',' '))
%     % 
%     % end
% 
% 
% end
% legend('No Empty - Test','Empty - Test')
% 
% % Let's manually increase ribosome binding affinity to see whether that
% % fixes the problem 
% 
%     % Create a simFunction for this 
% sfGFP_idx = find(strcmp(all_species_name_test,'protein sfGFP*')); 
% observables = [all_species_name_test(1:sfGFP_idx - 1),all_species_name_test(sfGFP_idx + 1:end)]; 
% simFunction_dataNames = [{'protein sfGFP*'},observables]; 
% simFunction = create_simFun_from_problemObject(RFM_integration_test_result_file.problemObject,observables); 
% dosing_information = create_dosing_info_from_problemObject(RFM_integration_test_result_file.problemObject); 
% 
%     % Modify ribosome binding affinity and solve simFunction 
% opt_fitResult_test = RFM_integration_test_result_file.all_fitResults{sorted_SSE_idx_test(1)}; 
% estimated_params = [opt_fitResult_test.ParameterEstimates.Estimate];
% estimated_param_names = {opt_fitResult_test.ParameterEstimates.Name};
% modify_param_name_list = {'TXTL_UTR_GFP_R','TXTL_UTR_kanR_R','TXTL_TL_init_k'}; 
% modify_target_value_list = [1e-09,1e-04,1e+04]; 
% for modify_idx = 1:length(modify_param_name_list)
%     modify_param_name = modify_param_name_list{modify_idx}; 
%     param_idx = find(strcmp(modify_param_name,estimated_param_names));
%     estimated_params(param_idx) = modify_target_value_list(modify_idx); 
% end
% [simulated_time,simulated_data] = simFunction(estimated_params',tEnd,dosing_information,tStart:tEnd);

%     % Plot out species time course - compare with unmodified parameters 
% figure; 
% for species_idx = 1:length(simFunction_dataNames)
% 
%     subplot(ceil(sqrt(length(simData_test_no_empty_conc_vec1.DataNames))),ceil(sqrt(length(simData_test_no_empty_conc_vec1.DataNames))),species_idx)
% 
%     species_name = simFunction_dataNames{species_idx};
%     simData_idx = find(strcmp(species_name,simData_test_no_empty_conc_vec1.DataNames)); 
% 
%     modified_param_no_empty_time = simulated_time{1}; 
%     modified_param_no_empty_data = simulated_data{1};
%     modified_param_empty_time = simulated_time{4};
%     modified_param_empty_data = simulated_data{4}; 
% 
%     plot(simData_test_no_empty_conc_vec1.Time,simData_test_no_empty_conc_vec1.Data(:,simData_idx),'LineWidth',1.5,'Color','g','LineStyle','--')
%     hold on 
%     plot(simData_test_empty_conc_vec1.Time,simData_test_empty_conc_vec1.Data(:,simData_idx),'LineWidth',1.5,'Color',[0.5,0.5,0.5],'LineStyle','--')
%     plot(modified_param_no_empty_time,modified_param_no_empty_data(:,species_idx),'LineWidth',1.5,'Color','g')
%     plot(modified_param_empty_time,modified_param_empty_data(:,species_idx),'LineWidth',1.5,'Color',[0.5,0.5,0.5])
% 
%     title(strrep(species_name,'_',' '))
% 
% end
% legend('RFM no empty','RFM empty','RFM no empty modified params','RFM empty modified params')
%% Troubleshoot sampled parameters in run 1000 & negative values in observed species 
% param_and_init_cond_sampling_result = load('test_save_files/20240613_param_init_cond_sampling_constrained.mat');
% satisfied_sampled_params = param_and_init_cond_sampling_result.nominal_sampled_params(param_and_init_cond_sampling_result.final_satisfied_idx_list,:); 

%     % Let's also manually check the parameter constraints again 
% param_names = {param_and_init_cond_sampling_result.problemObject.Estimated.Name}; 
% 
%     % Check (1) RNA deg kc (2) RNase binding (3) RNAP binding 
% RNA_deg_kc_related_params = {'TXTL_RNAdeg_kc_sfGFP','TXTL_RNAdeg_kc_kanR','TXTL_RNAdeg_kc'}; 
% RNase_binding_related_params = {'TXTL_RNAdeg_R_sfGFP','TXTL_RNAdeg_R'};
% RNAP_binding_related_params = {'TXTL_PT7_RNAPbound_R','TXTL_PT773_RNAPbound_R','TXTL_PJ23119_RNAPbound_R','TXTL_PJ23105_RNAPbound_R'};
% [RNA_deg_kc_param_vals,~] = extract_sampled_params(param_names,satisfied_sampled_params,RNA_deg_kc_related_params); 
% [RNase_binding_param_vals,~] = extract_sampled_params(param_names,satisfied_sampled_params,RNase_binding_related_params); 
% [RNAP_binding_param_vals,~] = extract_sampled_params(param_names,satisfied_sampled_params,RNAP_binding_related_params); 
% 
%     % RNA deg kc 
% RNA_deg_kc_sanity_check_1 = all(RNA_deg_kc_param_vals(:,1) < RNA_deg_kc_param_vals(:,3)); 
% RNA_deg_kc_sanity_check_2 = all(RNA_deg_kc_param_vals(:,2) < RNA_deg_kparam_vals(:,3)); 
%     % RNase binding 
% RNase_binding_sanity_check_1 = all(RNase_binding_param_vals(:,1) > RNase_binding_param_vals(:,2));
%     % RNAP binding 
% RNAP_binding_sanity_check_1 = all(RNAP_binding_param_vals(:,1) < RNAP_binding_param_vals(:,2)); 
% RNAP_binding_sanity_check_3 = all(RNAP_binding_param_vals(:,1) < RNAP_binding_param_vals(:,4)); 
% RNAP_binding_sanity_check_4 = all(RNAP_binding_param_vals(:,3) < RNAP_binding_param_vals(:,4)); 
% 
% function [sampled_param_vals,sampled_param_names] = extract_sampled_params(all_param_names,all_sampled_params,param_name_oi_list)
%     % Extract parameter values from param_name_oi_list
%         % Preassign 
%     sampled_param_vals = nan(size(all_sampled_params,1),length(param_name_oi_list)); 
%     sampled_param_names = cell(length(param_name_oi_list),1); 
% 
%     for param_idx = 1:length(param_name_oi_list)
%         param_name_oi = param_name_oi_list{param_idx}; 
%         param_name_oi_idx = strcmp(param_name_oi,all_param_names); 
%         sampled_param_vals(:,param_idx) = all_sampled_params(:,param_name_oi_idx); 
%         sampled_param_names{param_idx} = param_name_oi; 
%     end
% 
% end

%     % Find cases where residual mRNA conc penalty is smaller than 0 
% result_file = load('test_save_files/20240619_large_scale_parameter_sampling_summary.mat'); 
% all_sampled_params = result_file.all_sampled_params; 
% all_penalty_terms = result_file.all_penalty_terms; 
% all_dosing_mat = result_file.all_dosing_mat; 
% penalty_term_labels = result_file.penalty_term_labels; 
% residue_mRNA_penalty_term = all_penalty_terms(:,length(penalty_term_labels)); 
% negative_idx_list = find(residue_mRNA_penalty_term < 0); 
% negative_mRNA_penalty_term = residue_mRNA_penalty_term(negative_idx_list); 
% 
%     % Get sampled params & init cond to calculate the time-course
% sampled_params_oi_log = all_sampled_params(negative_idx_list,:);
% sampled_params_oi = 10 .^ sampled_params_oi_log; 
% dosing_information_oi_idx = mod(negative_idx_list,length(result_file.all_dosing_information)); 
% dosing_information_oi = result_file.all_dosing_information{dosing_information_oi_idx}; 
% 
%     % Create SimFunction first
%         % Calculate time-course data using the current set of parameters 
% dosing_target_list = param_and_init_cond_sampling_result.target_name_list;
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject = setProblemObject_v2(group_description,param_info_path,[]);
%         % Modify the track species list (making GFP the first species)
% model_species_name_list = {problemObject.Model.Species.Name}; 
% sfGFP_idx = find(strcmp(model_species_name_list,'protein sfGFP*')); 
% simFunction_tracked_species_list = [model_species_name_list(sfGFP_idx),model_species_name_list(1:sfGFP_idx - 1),model_species_name_list(sfGFP_idx + 1:end)]; 
% simFunction = createSimFunction(problemObject.Model,param_names,simFunction_tracked_species_list,dosing_target_list);
% 
%     % Solve the simulated data
% [simulated_time,simulated_data] = simFunction(sampled_params_oi,tEnd,dosing_information_oi,tStart:tEnd);
% % plot_simulated(simulated_time,simulated_data,'baseline')
% % plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% simulated_time_oi = simulated_time{3,1}; 
% simulated_data_oi = simulated_data{3,1}; 
% [rxn_name_list,rxn_time_course] = get_reaction_time_course(problemObject.Model,...
%     {problemObject.Estimated.Name},sampled_params_oi,simulated_time_oi,simulated_data_oi,simFunction_tracked_species_list); 
% figure; 
% for species_idx = 1:length(simFunction_tracked_species_list)
%     tracked_species_name = simFunction_tracked_species_list{species_idx}; 
%     tracked_species_for_plot = strrep(tracked_species_name,'_',' '); 
%     subplot(10,10,species_idx)
%     plot(simulated_time_oi,simulated_data_oi(:,species_idx))
%     title(tracked_species_for_plot); 
% end

%% Perform analysis on parameter sampling results 

%     % Load summary of sampling results 
% result_file = load('test_save_files/20240619_large_scale_parameter_sampling_summary.mat'); 
% all_sampled_params = result_file.all_sampled_params; 
% all_penalty_terms = result_file.all_penalty_terms; 
% all_dosing_mat = result_file.all_dosing_mat; 
% 
%     % Filter out cases where large positive crosstalk ratio is satisfied 
% large_positive_crosstalk_ratio_penalty = all_penalty_terms(:,5); 
% satisfied_idx_list = find(large_positive_crosstalk_ratio_penalty < -2); % note that the threshold is arbitrarily set to -1 for now  
% 
%     % In these cases, what other penalties are violated 
% baseline_data_penalty = all_penalty_terms(satisfied_idx_list,3); 
% transition_penalty = all_penalty_terms(satisfied_idx_list,4); 
% residue_mRNA_penalty = all_penalty_terms(satisfied_idx_list,6);
% 
% num_baseline_violated = sum(baseline_data_penalty > 0); 
% num_transition_violated = sum(transition_penalty > 0); 
% num_residue_mRNA_violated = sum(residue_mRNA_penalty > 10000); 
% 
% perc_violated = [num_baseline_violated,num_transition_violated,num_residue_mRNA_violated] ./ length(satisfied_idx_list); 
% 
% Xlabel = categorical({'Baseline','Crosstalk Transition','Residue mRNA'}); 
% Xlabel = reordercats(Xlabel,{'Baseline','Crosstalk Transition','Residue mRNA'});
% bar(Xlabel,perc_violated)

%     % Load the initially sampled parameters and initial conditions
% param_and_init_cond_sampling_result = load('test_save_files/20240613_param_init_cond_sampling_constrained.mat');
% 
%     % Load an example result file to get parameter names 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject = setProblemObject_v2(group_description,param_info_path,[]);
% param_names = {problemObject.Estimated.Name}; 
% 
% % Preprocess to remove missing values if necessary? 
% invalid_row_idx_list = any(isnan(all_penalty_terms),2); 
% valid_all_penalty_terms = all_penalty_terms(~invalid_row_idx_list,:); 
% valid_all_sampled_params = all_sampled_params(~invalid_row_idx_list,:); 



% % 1.1 distribution of penalty terms first 
% for penalty_idx = size(all_penalty_terms,2):size(all_penalty_terms,2)
%     penalty_term_label = result_file.penalty_term_labels{penalty_idx};  
% 
%     if mean(all_penalty_terms(:,penalty_idx)) > 10
%         penalty_term_for_plot = log10(all_penalty_terms(:,penalty_idx));
%         xlabel_name = sprintf('log10(%s)',strrep(penalty_term_label,'_',' ')); 
%     else
%         penalty_term_for_plot = all_penalty_terms(:,penalty_idx);
%         xlabel_name = strrep(penalty_term_label,'_',' '); 
%     end
%     figure; 
%     histogram(log10(penalty_term_for_plot(penalty_term_for_plot > 0)))
%     % histogram(penalty_term_for_plot(penalty_term_for_plot < 0))
%     xlabel(xlabel_name)
%     ylabel('# Occurence')
%     title(strrep(penalty_term_label,'_',' ')); 
%     set(gca,'FontSize',12); 
% 
% end

% % % 1.2 parameter distribution for satisfied vs.
% % unsatisfeid for each penalty
% satisfied_penalty_cut_off_list = [nan,nan,0,0,-1,100]; 
% 
% figure; 
% for param_idx = 1:size(all_sampled_params,2)
%     param_name = param_names{param_idx}; 
% 
%     subplot(ceil(sqrt(size(all_sampled_params,2))),ceil(sqrt(size(all_sampled_params,2))),param_idx)
% 
% 
%     for penalty_idx = 1:size(all_penalty_terms,2)
% 
%         % Get penalty term labels for plot titles & cut-off values 
%         penalty_term_label = result_file.penalty_term_labels{penalty_idx}; 
%         satisfied_penalty_cut_off = satisfied_penalty_cut_off_list(penalty_idx); 
% 
%         % Separate sampled parameters into satisfied vs. unsatisfied group 
%         satisfied_sampled_penalty_term_idx = valid_all_penalty_terms(:,penalty_idx) < satisfied_penalty_cut_off;
%         unsatisfied_sampled_penalty_term_idx = valid_all_penalty_terms(:,penalty_idx) > satisfied_penalty_cut_off;
% 
%         satisfied_sampled_params = valid_all_sampled_params(satisfied_sampled_penalty_term_idx,param_idx);
%         unsatisfied_sampled_params = valid_all_sampled_params(unsatisfied_sampled_penalty_term_idx,param_idx);
% 
%         % Plot out distribution in log space
%         histogram(satisfied_sampled_params)
%         hold on 
%         histogram(unsatisfied_sampled_params)
%         title(strrep(param_name,'_',' '))
% 
%     end
% 
% end 
% legend('Satisfied','Unsatified')

% % 2.2 Correlation-based sensitivity analysis. using PCC & PRCC - with parameters and initial conditions
% % separately 
%     % Note that these are recommended when the relationship is monotonic,
%     % which may or may not be true here 
% 
% PRCC = partialcorri(valid_all_penalty_terms,valid_all_sampled_params); 
% figure; 
% h1 = heatmap(strrep(param_names,'_',' '),strrep(result_file.penalty_term_labels,'_',' '),PRCC,'Colormap',jet); 
% xlabel('Parameter Names')
% ylabel('Penalty Terms') 
% 
% % Let's try using PCC & PRCC - with both parameters and initial conditions 
% mega_dosing_mat = repmat(param_and_init_cond_sampling_result.nominal_sampled_empty_init_cond,10000,1); 
% valid_mega_dosing_mat = mega_dosing_mat(~invalid_row_idx_list,:);
% augmented_param_initCond_mat = [valid_all_sampled_params,valid_mega_dosing_mat];
% augmented_param_initCond_labels = [strrep(param_names,'_',' '),{'Empty Plasmid Concentration'}]; % Param names and initial cond names 
% augmented_PCC = partialcorri(valid_all_penalty_terms,augmented_param_initCond_mat); 
% rank_valid_all_penalty_terms = tiedrank(valid_all_penalty_terms);
% rank_augmented_param_initCond_mat = tiedrank(augmented_param_initCond_mat); 
% augmented_PRCC = partialcorri(rank_valid_all_penalty_terms,rank_augmented_param_initCond_mat); 
% 
% figure; 
% h2 = heatmap(augmented_param_initCond_labels,strrep(result_file.penalty_term_labels,'_',' '),augmented_PCC,'Colormap',jet); 
% xlabel('Parameter & Initial Condition Names')
% ylabel('Penalty Terms') 
% title('Partial Correlation Coefficient')
% 
% figure; 
% h3 = heatmap(augmented_param_initCond_labels,strrep(result_file.penalty_term_labels,'_',' '),augmented_PRCC,'Colormap',jet); 
% xlabel('Parameter & Initial Condition Names')
% ylabel('Penalty Terms') 
% title('Partial Rank Correlation Coefficient')

% % What about calculating just correlation? 
% corr_matrix = corr(valid_all_penalty_terms,augmented_param_initCond_mat); 
% figure;
% h4 = heatmap(augmented_param_initCond_labels,strrep(result_file.penalty_term_labels,'_',' '),corr_matrix,'Colormap',jet); 
% xlabel('Parameter & Initial Condition Names')
% ylabel('Penalty Terms') 

% Is there a sample that satisfies all phenotype criteria? 

%     % We don't really care about deviations 
% satisfied_idx_list = valid_all_penalty_terms(:,3) < 0 & valid_all_penalty_terms(:,4) < 0 & valid_all_penalty_terms(:,5) < -1 & valid_all_penalty_terms(:,6) < 100; 
%     % Check the sampled parameter values and initial condition 
% satisfied_sampled_params = valid_all_sampled_params(satisfied_idx_list,:); 
% satisfied_init_cond = valid_mega_dosing_mat(satisfied_idx_list,:); 
% nominal_satisfied_sampled_params = 10 .^ satisfied_sampled_params;
% nominal_satisfied_init_cond = 10 .^ satisfied_init_cond; 

%     % Empty Plasmid Concentration Distribution Visulization
% figure; 
% histogram(log10(satisfied_init_cond),'FaceColor',"#D95319")
% title('Empty Plasmid Concentration Satisfying Qualitative Phenotypes')
% xlabel('log_{10}([EmptyPlasmid])')
% ylabel('# Occurence')
% set(gca,'FontSize',12)
% 
%     % Parameter Distribution Visualization
% figure;
% for param_idx = 1:length(param_names)
%     selected_satisfied_sampled_params = satisfied_sampled_params(:,param_idx); 
%     param_name = param_names{param_idx}; 
%     subplot(ceil(sqrt(length(param_names))),ceil(sqrt(length(param_names))),param_idx)
%     histogram(selected_satisfied_sampled_params)
%     title(strrep(param_name,'_',' '));
%     xlabel('log_{10}parameter')
% end
% sgtitle('Parameter Distribution Satisfying Qualitative Phenotypes');

    % What does the deviation look like? 
% satisfied_penalty_terms = valid_all_penalty_terms(satisfied_idx_list,:);

%     % Plot out 10 of these 
%         % Create SimFunction first
%     % Calculate time-course data using the current set of parameters 
% dosing_target_list = param_and_init_cond_sampling_result.target_name_list;
%     % Modify the track species list (making GFP the first species)
% model_species_name_list = {problemObject.Model.Species.Name}; 
% sfGFP_idx = find(strcmp(model_species_name_list,'protein sfGFP*')); 
% simFunction_tracked_species_list = [model_species_name_list(sfGFP_idx),model_species_name_list(1:sfGFP_idx - 1),model_species_name_list(sfGFP_idx + 1:end)]; 
% simFunction = createSimFunction(problemObject.Model,param_names,simFunction_tracked_species_list,dosing_target_list); 
% 
% num_satisfied_samples = sum(satisfied_idx_list); 
% temp_sequence = randperm(num_satisfied_samples); 
% for i = 1:10
%     selected_idx = temp_sequence(i); 
%     sampled_params_oi = nominal_satisfied_sampled_params(selected_idx,:); 
%     sampled_penalties_oi = satisfied_penalty_terms(selected_idx,:); 
%     sampled_init_cond_oi = satisfied_init_cond(selected_idx,:);
%     init_cond_oi_idx = find(param_and_init_cond_sampling_result.nominal_sampled_empty_init_cond == sampled_init_cond_oi); 
%     dosing_info_oi = result_file.all_dosing_information{init_cond_oi_idx}; 
% 
%     [simulated_time,simulated_data] = simFunction(sampled_params_oi,tEnd,dosing_info_oi,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% 
% end


%% Organize parameter sampling results 
% num_params_iter = 1000;
% num_dosing = 100; 
% num_total_iter = num_params_iter * num_dosing; 
% 
% for param_iter = 1:num_params_iter
% 
%     result_file_name = sprintf('20240613_param_and_init_cond_sampling_constrained_run%d.mat',param_iter);
%     try
%         result_file = load(result_file_name); 
%             % Initialize when first loading a result file successfully 
%         if ~exist('penalty_term_labels','var') 
%             penalty_term_labels = result_file.penalty_labels; 
%             all_dosing_information = result_file.all_dosing_information; 
%             all_dosing_mat = result_file.all_dosing_mat; 
%             num_params = size(result_file.sampled_params_selected,2); 
%             all_sampled_params = nan(num_total_iter,num_params); 
%             all_penalty_terms = nan(num_total_iter,length(penalty_term_labels)); 
%         end
% 
%             % Record penalty terms
%         for local_param_iter = 1:size(result_file.sampled_params_selected,1)
%             start_idx = (param_iter - 1) * size(result_file.sampled_params_selected,1) * num_dosing + (local_param_iter - 1) * num_dosing + 1; 
%             stop_idx = (param_iter - 1) * size(result_file.sampled_params_selected,1) * num_dosing + local_param_iter * num_dosing;
%             sampled_params = result_file.sampled_params_selected(local_param_iter,:);
%             all_sampled_params(start_idx:stop_idx,:) = repmat(sampled_params,size(all_dosing_mat,1),1);
%             all_penalty_terms(start_idx:stop_idx,:) = result_file.all_penalty_terms((local_param_iter - 1) * num_dosing + 1:local_param_iter * num_dosing,:);
%         end
% 
% 
%     catch
%         fprintf('\n iter #%d does not exist',param_iter); 
%     end
% 
% 
% end

% save('test_save_files/20240619_large_scale_parameter_sampling_summary.mat');

%% TX gen TL toxin mechanism failure mode check 
% % result_file_name = 'param_est_run_save/20230125_param_est_run1111_TXgen_Ribo.mat'; 
% % result_file = load(result_file_name); 
% all_fitResults = result_file.all_fitResults; 
% % figure;
% % for iter = 1:length(all_fitResults)
% %     fitResult_oi = all_fitResults{iter}; 
% %     fitted_data = fitted(fitResult_oi); 
% % 
% %     [Time,Data] = process_crosstalk_data_from_source(fitted_data,'simulated_data');
% %     % [Time,Data] = process_crosstalk_data_from_source(problemObject.Data,'grouped_data');
% %     subplot(7,7,iter)
% %     for conc_idx = 1:length(Data)
% %         plot(Time{conc_idx}./3600,Data{conc_idx},'LineWidth',1.5)
% %         hold on 
% %     end
% %     sgtitle('TX-gen toxin on Ribo')
% %     % xlabel('Time (hr)','FontSize',14)
% %     % ylabel('Concentration (nM)','FontSize',14)
% %     % legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')
% %     set(gca,'FontSize',14)
% % end
% 
% % Select iter #33 as a case study 
% fitResult_oi = all_fitResults{33,1}; 
% 
%     % Create simFunction and track all species
% problemObject = result_file.problemObject; 
% all_species_name_list = {problemObject.Model.Species.Name};
% GFP_species_idx = find(strcmp({problemObject.Model.Species.Name},'protein sfGFP*'));
% additional_track_species = [all_species_name_list{1:GFP_species_idx-1},all_species_name_list(GFP_species_idx+1:end)]; 
% species_name_list = [{'GFP'},additional_track_species];
% simFunction = create_simFun_from_problemObject(problemObject,additional_track_species); 
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
%     % Solve the simFunction 
% estimated_params = [fitResult_oi.ParameterEstimates.Estimate]; 
%     % Modify estimated params 
% estimated_params(strcmp({problemObject.Estimated.Name},'TXTL_RNAdeg_R_sfGFP')) = 1;  
% estimated_params(strcmp({problemObject.Estimated.Name},'k_toxin')) = 0.001; 
% [simulated_time,simulated_data] = simFunction(estimated_params',tEnd,dosing_information,tStart:tEnd);    
% 
% figure;
% for conc_idx = 1:length(simulated_data)
%     plot(simulated_time{conc_idx}./3600,simulated_data{conc_idx}(:,1),'LineWidth',1.5)
%     hold on 
% end
% sgtitle('TX-gen toxin on Ribo')
% xlabel('Time (hr)','FontSize',14)
% ylabel('Concentration (nM)','FontSize',14)
% legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')
% set(gca,'FontSize',14)
% 
%     % Organize solution into a compatible format 
% simData.Time = simulated_time{7}; 
% simData.Data = simulated_data{7}; 
% simData.DataNames = [{'protein sfGFP*'},additional_track_species]; 
% [rxn_name_list,rxn_time_course] = get_reaction_time_course(problemObject.Model,{problemObject.Estimated.Name},estimated_params,...
%     simData.Time,simData.Data,simData.DataNames);
% 
%     % Plot out time-course for key species 
% reporter_RNA_species_idx_list = contains(simData.DataNames,'RNA') & ~contains(simData.DataNames,'RNase') & contains(simData.DataNames,'sfGFP'); 
% reporter_RNA_species = simData.DataNames(reporter_RNA_species_idx_list); 
% kanR_RNA_species_idx_list = contains(simData.DataNames,'RNA') & ~contains(simData.DataNames,'RNase') & contains(simData.DataNames,'kanR'); 
% kanR_RNA_species = simData.DataNames(kanR_RNA_species_idx_list);
% species_oi_list = [{'AGTP','CUTP','RNAP','t7RNAP','Ribo','RNase','protein kanR','toxin'},reporter_RNA_species,kanR_RNA_species]; 
% species_oi_name_list = {'AGTP','CUTP','RNAP','t7RNAP','Ribo','RNase','protein kanR','toxin','sfGFP RNA','kanR RNA'}; 
% species_oi_idx_list = [1:9,length(species_oi_list)-length(kanR_RNA_species) + 1];
% [timeVec,species_time_course] = get_species_time_course(simData,species_oi_list);
% figure; 
% for species_oi_idx = 1:length(species_oi_name_list)
%     species_oi_name = species_oi_name_list{species_oi_idx};
%     species_oi_actual_start_idx = species_oi_idx_list(species_oi_idx); 
%     if ~isequal(species_oi_idx,length(species_oi_name_list))
%         species_oi_actual_end_idx = species_oi_idx_list(species_oi_idx + 1); 
%     else
%         species_oi_actual_end_idx = length(species_oi_name_list); 
%     end
%     species_oi_time_course = sum(species_time_course(:,species_oi_actual_start_idx:species_oi_actual_end_idx),2);
%     subplot(ceil(sqrt(length(species_oi_name_list))),ceil(sqrt(length(species_oi_name_list))),species_oi_idx)
%     plot(timeVec,species_oi_time_course)
%     xlabel('Time(s)')
%     title(species_oi_name_list{species_oi_idx}) 
% end

%% Initial sampling analysis
% meta_all_penalty_terms = nan(1e+8,6); 
% meta_all_sampled_params = nan(10000,31); 
% for run_idx = 1:150
%     result_file_name = sprintf('param_est_run_save/20240612_param_and_init_cond_sampling_run%d.mat',run_idx); 
%     try
%         result_file = load(result_file_name); 
%         all_penalty_terms = result_file.all_penalty_terms; 
%         all_sampled_params = result_file.sampled_params_selected; 
%         all_dosing_mat = result_file.all_dosing_mat; 
%         meta_all_penalty_terms((run_idx - 1) * 1000 + 1:run_idx * 1000,:) = all_penalty_terms; 
%         meta_all_sampled_params((run_idx - 1) * 10 + 1:run_idx * 10,:) = all_sampled_params; 
%     catch
%         fprintf('\n Sampling run %d not done',run_idx)
%     end
% 
% end

% 
% penalty_term_labels = result_file.penalty_labels; 
% % Plot penalty distributions 
% for penalty_idx = 1:length(penalty_term_labels)
%     penalty_term_label = penalty_term_labels{penalty_idx}; 
%     penalty_term_oi = meta_all_penalty_terms(:,penalty_idx); 
%     penalty_term_oi_for_plot = penalty_term_oi(~isnan(penalty_term_oi)); 
%     figure;
%     histogram(penalty_term_oi);
%     title(penalty_term_label)
%     ylabel('#Occurence')
% end

% Filter for rows with all constraints satisfied 
% satisfied_idx_list = meta_all_penalty_terms(:,3) < 0 & meta_all_penalty_terms(:,4) < 0 & meta_all_penalty_terms(:,5) < 0 & log10(meta_all_penalty_terms) < 2; 


%% Filter out sampled parameters based on existing parameter constraints 

% load('test_save_files/20240612_param_init_cond_sampling.mat'); 
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


%% Parameter + initial condition sampling 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject = setProblemObject_v2(group_description,param_info_path,[]); 
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
% % Fix mRNA degradation parameters (using 95% CI as LB/UB) 
% mRNA_deg_param_info_file = load('test_save_files/20240530_mRNA_degradation_param_bounds.mat'); 
% for estimated_param_idx = 1:length(mRNA_deg_param_info_file.estimated_params_labels)
%     estimated_param_label = mRNA_deg_param_info_file.estimated_params_labels{estimated_param_idx};
%     estimated_param_bounds = mRNA_deg_param_info_file.param_bounds_CI(estimated_param_idx,:);
% 
%     kinetic_param_idx = strcmp(estimated_param_label,kinetic_param_names); 
%     kinetic_param_bounds(kinetic_param_idx,:) = estimated_param_bounds; 
% 
% end
% log_param_bounds = log10(kinetic_param_bounds); 
% 
% % Generate sampled kinetic parameter values 
% num_sample = 400000; 
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
% 
% % Generate sampled plasmid concentrations 
% num_init_cond_sample = 100; 
%     % Sample initial condition for empty plasmid first 
% log_sampled_empty_init_cond = -2 + 4 * lhsdesign(num_init_cond_sample,1); 
% nominal_sampled_empty_init_cond = 10 .^ log_sampled_empty_init_cond; 
% %     % Sample initial conditions for T7 & sigma70 reporter plasmids 
% % lhsSamples_init_cond = lhsdesign(num_init_cond_sample,1); 
% % sampled_init_cond = -2 + 4 * lhsSamples_init_cond; % Sample initial plasmid concentration from the [0.01,100] interval 
% % nominal_sampled_init_cond = 10 .^ sampled_init_cond; 
% 
% % Create dosing information based on sampled initial condition
% 
%     % Get target names of interest for PE 
% model_species_list = {modelObj.Species.Name}; 
% T7_reporter_target_idx_list = contains(model_species_list,'DNA') & contains(model_species_list,'T7') & contains(model_species_list,'sfGFP') & ~contains(model_species_list,':'); 
% T7_reporter_target_list = model_species_list(T7_reporter_target_idx_list);
% sigma70_reporter_target_idx_list = contains(model_species_list,'DNA') & contains(model_species_list,'J23') & contains(model_species_list,'sfGFP') & ~contains(model_species_list,':'); 
% sigma70_reporter_target_list = model_species_list(sigma70_reporter_target_idx_list);
% empty_target_idx_list = contains(model_species_list,'DNA') & (contains(model_species_list,'utrempty') | contains(model_species_list,'utrkanR')) & ~contains(model_species_list,':'); 
% empty_target_list = model_species_list(empty_target_idx_list);
% target_name_list = [T7_reporter_target_list,sigma70_reporter_target_list,empty_target_list]; 
% 
%     % Pre-code the dosing information object 
% all_dosing_mat = cell(num_init_cond_sample,1); 
% all_dosing_information = cell(num_init_cond_sample,1); 
% 
%     % Convert sampled initial conditions into a dosing matrix 
%         % Each dosing matrix - 1 initial empty * sampled reporter conc * 2 
% for init_cond_idx_1 = 1:num_init_cond_sample
%     nominal_sampled_empty_init_cond_single = nominal_sampled_empty_init_cond(init_cond_idx_1); 
% 
%     % Assign initial conditions to dosing information
%     dosing_information_mat = zeros(length(conc_vec) * 16,length(target_name_list)); 
%     for target_idx = 1:length(T7_reporter_target_list) + length(sigma70_reporter_target_list)
% 
%         for comb_idx = 1:num_comb
% 
%             start_idx = (target_idx - 1) * num_comb * length(conc_vec) + (comb_idx - 1) * length(conc_vec) + 1; 
%             end_idx =  (target_idx - 1) * num_comb * length(conc_vec) + comb_idx * length(conc_vec); 
%             dosing_information_mat(start_idx:end_idx,target_idx) = conc_vec;
%             if ~isequal(comb_idx,1)
%                     dosing_information_mat(start_idx:end_idx,length(T7_reporter_target_list) + length(sigma70_reporter_target_list) + comb_idx - 1) = ...
%                 nominal_sampled_empty_init_cond_single;
%             end
%         end
% 
%     end
% 
%     % Convert dosing matrix into a Doses object 
%     dosing_information = cell(size(dosing_information_mat)); 
%     for row_idx = 1:size(dosing_information_mat,1)
%         for col_idx = 1:size(dosing_information_mat,2)
%             dosing_table = array2table([0,dosing_information_mat(row_idx,col_idx)]); 
%             dosing_table.Properties.VariableNames = {'Time','Amount'};
%             dosing_information{row_idx,col_idx} = dosing_table; 
%         end
%     end
% 
%     all_dosing_mat{init_cond_idx_1,1} = dosing_information_mat; 
%     all_dosing_information{init_cond_idx_1,1} = dosing_information;
% 
% end

% % Evaluate model and calculate objective values 
% num_penalty_terms = 6; 
% all_penalty_terms = nan(num_sample * num_init_cond_sample,num_penalty_terms); 
% 
% for sample_idx = 1:num_sample
%     sampled_params_single = sampled_params(sample_idx,:); 
%     nominal_sampled_params_single = 10 .^ sampled_params_single; 
% 
%     % Sample initial conditions and create new dosing information 
%     for init_sample_idx = 1:num_init_cond_sample
%         dosing_information = all_dosing_information{init_sample_idx}; 
% 
%         % Evaluate penalty values 
%         [penalty_values,penalty_labels] = wrapper_calculate_obj_v2(modelObj,experimental_grouped_data,nominal_sampled_params_single,kinetic_param_names,dosing_information,target_name_list); 
% 
%         all_penalty_terms((sample_idx - 1) * num_init_cond_sample + init_sample_idx,:) = penalty_values; 
%     end
% 
% end


%% Analyze parameter distribution and SSE distribution with large-scale fitting 
% large_scale_fitting_result_file = load('param_est_run_save/20231106_large_scale_run.mat');
% large_scale_fitting_result_struct = large_scale_fitting_result_file.updated_optimization_run_table; 
% 
% %     % SSE distribution 
% % all_SSE = [large_scale_fitting_result_struct.SSE];
% % histogram(log10(all_SSE));
% 
%     % Parameter distribution 
% 
%         % Get basic information on the result struct 
% sample_estimated_params = cell2mat(large_scale_fitting_result_struct{1,'ParameterEstimates'});
% num_params = length(sample_estimated_params); 
% num_runs = height(large_scale_fitting_result_struct); 
% 
%         % Convert estimated params into the correct size 
% all_estimated_params_cell = [large_scale_fitting_result_struct.ParameterEstimates];
% all_estimated_params_vec = cell2mat(all_estimated_params_cell);
% all_estimated_params_mat = reshape(all_estimated_params_vec,[num_params,num_runs]);
% all_estimated_params_mat = all_estimated_params_mat'; 
% 
%     % Create a problemObject & get parameter names * bounds 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_test_PE.xlsx';
% problemObject = setProblemObject_PE(group_description,param_info_path,[]); 
% Estimated = problemObject.Estimated; 
% estimated_param_names = {Estimated.Name}; 
% param_bounds = {Estimated.Bounds}; 
% param_bounds_mat = cell2mat(param_bounds); 
% param_bounds_mat = reshape(param_bounds_mat,[2,num_params]); 
% param_bounds_mat = param_bounds_mat'; 
% 
% %     % Plot out distribution and bounds for each parameter 
% % figure; 
% % for param_idx = 1:num_params
% %     subplot(ceil(sqrt(num_params)),ceil(sqrt(num_params)),param_idx)
% %     % log transform each parameter? 
% %     estimated_params_vec = all_estimated_params_mat(:,param_idx); 
% %     param_bounds_single = param_bounds_mat(param_idx,:); 
% %     estimated_param_name = estimated_param_names{param_idx}; 
% %     estimated_param_name_for_plot = strrep(estimated_param_name,'_',' ');
% %     estimated_param_name_for_plot = strrep(estimated_param_name_for_plot,'TXTL',''); 
% %     log_estimated_params_vec = log10(estimated_params_vec); 
% %     histogram(log_estimated_params_vec);
% %     xline(log10(param_bounds_single),'-',{'LB','UB'})
% %     title(estimated_param_name_for_plot)
% % 
% % end
% 
%     % Evaluate other penalty terms in the large-scale parameter fitting &
%     % add to table 
% 
%     % Define inputs for penalty terms calculation 
% obj_string = all_obj_string{12}; 
% experimental_grouped_data = problemObject.Data;
% species_name_list = {problemObject.Model.Species.Name}; 
% rxn_name_list = {problemObject.Model.Reactions.Reaction};
% remove_idx = strcmp('protein sfGFP*',species_name_list); 
% additional_track_species_list = species_name_list(~remove_idx); 
% simulated_data_name_list = [{'protein sfGFP*'},additional_track_species_list];
% 
% 
% rxn_oi_list = {'[CUTP:AGTP:RNAP:DNA pJ23105--utrGFP--sfGFP] -> [term_RNAP:DNA pJ23105--utrGFP--sfGFP] + [RNA utrGFP--sfGFP]'}; 
% param_oi_list = {'TXTL_PJ23105_RNAPbound_R'}; 
% modify_value_magnitude_list = 0:2:6; 
% modify_values = 10.^ modify_value_magnitude_list;
% modify_value_list = {modify_values}; 
% sensitivity_goal_list = {1}; 
% sensitivity_options = struct('rxn_oi_list',rxn_oi_list,'param_oi_list',param_oi_list,'modify_value_list',modify_value_list,...
%     'sensitivity_goal_list',sensitivity_goal_list); 
% % sensitivity_options = []; 
% simFunction = create_simFun_from_problemObject(problemObject,additional_track_species_list); 
% 
% % all_obj_value = nan(size(all_estimated_params_mat,1),1); 
% % all_baseline_data_dev = nan(size(all_estimated_params_mat,1),1); 
% % all_crosstalk_ratio_dev = nan(size(all_estimated_params_mat,1),1); 
% % all_log_baseline_data_penalty = nan(size(all_estimated_params_mat,1),1); 
% % all_crosstalk_ratio_transition_penalty = nan(size(all_estimated_params_mat,1),1); 
% % all_crosstalk_ratio_promotor_strength_penalty = nan(size(all_estimated_params_mat,1),1); 
% all_promotor_sensitivity_penalty = nan(size(all_estimated_params_mat,1),1); 
% % all_other_penalty = nan(size(all_estimated_params_mat,1),1); 
% 
% for row_idx = 1:size(all_estimated_params_mat,1)
%     estimated_params_row = all_estimated_params_mat(row_idx,:); 
% 
%     [obj,baseline_data_dev,crosstalk_ratio_dev,log_baseline_data_penalty,crosstalk_ratio_transition_penalty,...
%         crosstalk_ratio_promotor_strength_penalty,promotor_sensitivity_penalty,other_penalty] = ...
%         wrapper_calculate_obj(estimated_params_row,obj_string,simFunction,problemObject,experimental_grouped_data,...
%         simulated_data_name_list,sensitivity_options,mode,num_conc,num_promotor);
% 
%     % all_obj_value(row_idx) = obj; 
%     % all_baseline_data_dev(row_idx) = baseline_data_dev; 
%     % all_crosstalk_ratio_dev(row_idx) = crosstalk_ratio_dev; 
%     % all_log_baseline_data_penalty(row_idx) = log_baseline_data_penalty; 
%     % all_crosstalk_ratio_transition_penalty(row_idx) = crosstalk_ratio_transition_penalty; 
%     % all_crosstalk_ratio_promotor_strength_penalty(row_idx) = crosstalk_ratio_promotor_strength_penalty; 
%     all_promotor_sensitivity_penalty(row_idx) = promotor_sensitivity_penalty; 
%     % all_other_penalty(row_idx) = other_penalty; 
% 
% end
% % [large_scale_fitting_result_struct.all_obj_value] = all_obj_value;
% % [large_scale_fitting_result_struct.all_baseline_data_dev] = all_baseline_data_dev;
% % [large_scale_fitting_result_struct.all_crosstalk_ratio_dev] = all_crosstalk_ratio_dev;
% % [large_scale_fitting_result_struct.all_log_baseline_data_penalty] = all_log_baseline_data_penalty;
% % [large_scale_fitting_result_struct.all_crosstalk_ratio_transition_penalty] = all_crosstalk_ratio_transition_penalty;
% % [large_scale_fitting_result_struct.all_crosstalk_ratio_promotor_strength_penalty] = all_crosstalk_ratio_promotor_strength_penalty;
% [large_scale_fitting_result_struct.all_promotor_sensitivity_penalty] = all_promotor_sensitivity_penalty;
% % [large_scale_fitting_result_struct.all_other_penalty] = all_other_penalty;
% save('test_save_files/20240531_large_scale_fitting_penalty_term.mat'); 

%% Check results for different obj strings 
%     % Outputs: Baseline expression + crosstalk ratios + penalty term
%     % breakdown for all iterations (highlighting top 5?) 
% directory_path = 'param_est_run_save'; 
% directory_files = dir(directory_path); 
% file_names = {directory_files.name};
% 
% for obj_string_idx = 15:15
% 
%     obj_string = all_obj_string{obj_string_idx}; 
% 
%     resultFileName_idx = contains(file_names,sprintf('objString_%d',obj_string_idx)); 
%     resultFileName = file_names{resultFileName_idx};    
%     resultFile = load(resultFileName); 
% 
%     % Check baseline and crosstalk ratio 
%     problemObject = resultFile.problemObject;
%     dosing_information = create_dosing_info_from_problemObject(problemObject); 
%     simFunction = resultFile.simFunction; 
%     all_fitResults = resultFile.all_fitResults; 
% 
%     all_fval = nan(length(all_fitResults),1); 
%     for iter = 1:length(all_fitResults)
%         fitResult = all_fitResults{iter,1}; 
%         all_fval(iter,1) = fitResult.fval; 
%     end
% 
%     [sorted_fval,sort_idx] = sort(all_fval,'ascend'); 
% 
%     for idx = 1:5
% 
%         fitResult = all_fitResults{sort_idx(idx),1}; 
%         estimated_params = fitResult.estimated_params; 
%         [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
%         plot_simulated(simulated_time,simulated_data,'baseline')
%         plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
%         sgtitle(sprintf('%s \n fval = %.2f',strrep(obj_string,'_',' '),sorted_fval(idx)))
% 
%     end
%     % save(resultFileName,'all_fval','sorted_fval','sort_idx','-append')
% 
%     % % Track penalty value change 
%     %     % Break obj string down into segments for eval statements 
%     % separate_idx_list = strfind(obj_string,'+'); 
%     % penalty_term_str_list = cell(length(resultFile.penalty_term_label),1);
%     % for sep_idx = 1:length(separate_idx_list) + 1 
%     %     if isequal(sep_idx,1)
%     %         start_idx = 1; 
%     %     else
%     %         start_idx = end_idx; 
%     %     end
%     %     if isequal(sep_idx,length(separate_idx_list) + 1)
%     %         end_idx = length(obj_string); 
%     %     else
%     %         end_idx = separate_idx_list(sep_idx); 
%     %     end
%     %     term_oi = obj_string(start_idx:end_idx);
%     %     if contains(term_oi,'^')
%     %         term_oi = strrep(term_oi,'^','.^');
%     %     end
%     %     penalty_term_str_list{sep_idx} = strrep(term_oi,'+','');
%     % end
% 
%         % Use a stacked bar plot in each case 
%     % figure; 
%     % all_penalty_for_plot = nan(3 * size(resultFile.all_penalty_term_final,1),length(resultFile.penalty_term_label));
%     % X_labels_for_plot = cell(3 * size(resultFile.all_penalty_term_final,1),1); 
%     % for iter = 1:size(resultFile.all_penalty_term_final,1)
%     %     penalty_term_init = resultFile.all_penalty_term_init(iter,:);
%     %     penalty_term_final = resultFile.all_penalty_term_final(iter,:); 
%     %     all_penalty_for_plot((iter - 1) * 3 + 1,:) = penalty_term_init; 
%     %     all_penalty_for_plot((iter - 1) * 3 + 2,:) = penalty_term_final;
%     %     X_labels_for_plot{(iter - 1) * 3 + 1} = sprintf('#%d init',iter); 
%     %     X_labels_for_plot{(iter - 1) * 3 + 2} = sprintf('#%d final',iter);
%     %     X_labels_for_plot{(iter - 1) * 3 + 3} = num2str(iter);  
%     % end
%     % X_labels = categorical(X_labels_for_plot); 
%     % X_labels = reordercats(X_labels,X_labels_for_plot); 
%     % 
%     % 
%     % modified_penalties = all_penalty_for_plot; 
%     % for penalty_idx = 1:length(penalty_term_str_list)
%     %     penalty_term_label = resultFile.penalty_term_label{penalty_idx}; 
%     %     eval(sprintf('%s = all_penalty_for_plot(:,penalty_idx);',penalty_term_label)); 
%     %     eval(sprintf('modified_penalty = %s;',penalty_term_str_list{penalty_idx})); 
%     %     modified_penalties(:,penalty_idx) = modified_penalty; 
%     % end
%     % bar(X_labels,modified_penalties,'stacked')
%     % legend(strrep(resultFile.penalty_term_label,'_',' '))
%     % sgtitle(sprintf('objective string %d',obj_string_idx))
% 
%     % for iter = 1:size(resultFile.all_penalty_term_final,1)
%     %     penalty_term_init = resultFile.all_penalty_term_init(iter,:);
%     %     penalty_term_final = resultFile.all_penalty_term_final(iter,:); 
%     %     penalties = [penalty_term_init;penalty_term_final];
%     %     modified_penalties = penalties; 
%     %     for penalty_idx = 1:length(penalty_term_str_list)
%     %         penalty_term_label = resultFile.penalty_term_label{penalty_idx}; 
%     %         eval(sprintf('%s = penalties(:,penalty_idx);',penalty_term_label)); 
%     %         eval(sprintf('modified_penalty = %s;',penalty_term_str_list{penalty_idx})); 
%     %         modified_penalties(:,penalty_idx) = modified_penalty; 
%     %     end
%     %     subplot(7,7,iter)
%     %     bar(modified_penalties,'stacked')
%     % end
%     % legend(strrep(resultFile.penalty_term_label,'_',' '))
%     % sgtitle(sprintf('objective string %d',obj_string_idx))
% 
%         % Add another plot to show any tradeoff among penalty terms 
% 
%     % For each initial points, show how many penalty terms
%     % increase/decrease in the optimization 
%     % figure; 
%     % num_increased_penalty = nan(size(resultFile.all_penalty_term_final,1),1);
%     % num_decreased_penalty = nan(size(resultFile.all_penalty_term_final,1),1);
%     % for iter = 1:size(resultFile.all_penalty_term_final,1)
%     %     penalty_term_init = resultFile.all_penalty_term_init(iter,:);
%     %     penalty_term_final = resultFile.all_penalty_term_final(iter,:); 
%     %     penalty_diff = penalty_term_final - penalty_term_init; 
%     % 
%     %     num_increased_penalty(iter) = sum(penalty_diff > 0);
%     %     num_decreased_penalty(iter) = sum(penalty_diff < 0);
%     % end
%     % num_decreased_penalty_for_plot = -1 .* num_decreased_penalty; 
%     % bar([num_increased_penalty,num_decreased_penalty_for_plot],'stacked')
%     % xlabel('Iteration')
%     % ylabel('# penalty terms')
% 
% end


%% Resolve the naming chaos
% %     % When running 48 iterations on cluster separately, I forgot to add the
% %     % iteration number into the saving path, so we got the same objString_X
% %     % with different run start time but don't know which iteration it is 
% % 
% % Loop through all result files 
% directory_path = 'param_est_run_save'; 
% directory_files = dir(directory_path); 
% file_names = {directory_files.name};
% 
% % covered_obj_string_idx_list = []; 
% % covered_obj_string_init_penalty = {}; 
% % all_penalty_term_init = nan(48,7); 
% for name_idx = 1:length(file_names)
%     file_name = file_names{name_idx}; 
%     if contains(file_name,'20240516') && contains(file_name,'objString_')
%         resultFile = load(file_name); 
% 
%         % Get obj string 
%         str_idx = strfind(file_name,'objString_'); 
%         obj_string_idx_format_str = file_name(str_idx+length('objString_'):str_idx+length('objString_')+1); 
%         obj_string_idx = str2double(obj_string_idx_format_str); 
%         obj_string = all_obj_string{obj_string_idx}; 
% 
%         % % Compile initial penalty term values to correspond to iteration #
%         % if ~ismember(obj_string_idx,covered_obj_string_idx_list)
%         % 
%         %     % Load variables required for the calculate penalty function 
%         %     simFunction = resultFile.simFunction; 
%         %     problemObject = resultFile.problemObject;
%         %     experimental_grouped_data = problemObject.Data; 
%         %     species_name_list = {problemObject.Model.Species.Name};
%         % 
%         %     for iter = 1:48
%         %         % Initial params for the 48 iterations
%         %         rng(iter)
%         %         bounds = [problemObject.Estimated.Bounds];
%         %         bounds_reorg = reshape(bounds,[2,numel(bounds)/2]);
%         %         lb = bounds_reorg(1,:);
%         %         ub = bounds_reorg(2,:);
%         %         init_params = lb + (ub - lb) .* rand(size(lb)); 
%         % 
%         %         % Define sensitivity options 
%         %         rxn_oi_list = {'[CUTP:AGTP:RNAP:DNA pJ23105--utrGFP--sfGFP] -> [term_RNAP:DNA pJ23105--utrGFP--sfGFP] + [RNA utrGFP--sfGFP]'}; 
%         %         param_oi_list = {'TXTL_PJ23105_RNAPbound_R'}; 
%         %         modify_value_magnitude_list = -4:1:4; 
%         %         modify_values = 10.^ modify_value_magnitude_list;
%         %         modify_value_list = {modify_values}; 
%         %         sensitivity_goal_list = {1}; 
%         %         sensitivity_options = struct('rxn_oi_list',rxn_oi_list,'param_oi_list',param_oi_list,'modify_value_list',modify_value_list,...
%         %             'sensitivity_goal_list',sensitivity_goal_list); 
%         % 
%         %         % Caclulate Initial penalty term values 
%         %         [~,baseline_data_dev_init,crosstalk_ratio_dev_init,log_baseline_data_penalty_init,crosstalk_ratio_transition_penalty_init,crosstalk_ratio_promotor_strength_penalty_init,promotor_sensitivity_penalty_init,other_penalty_init] = ...
%         %         wrapper_calculate_obj(init_params,obj_string,simFunction,problemObject,experimental_grouped_data,species_name_list,sensitivity_options,mode,num_conc,num_promotor); 
%         %         all_penalty_term_init(iter,:) = [baseline_data_dev_init,crosstalk_ratio_dev_init,log_baseline_data_penalty_init,crosstalk_ratio_transition_penalty_init,crosstalk_ratio_promotor_strength_penalty_init,promotor_sensitivity_penalty_init,other_penalty_init]; 
%         %     end
%         %     covered_obj_string_idx_list = [covered_obj_string_idx_list obj_string_idx];
%         %     covered_obj_string_init_penalty{end + 1} = all_penalty_term_init; 
%         % else
% 
%         % Identify the corresponding init penalty terms 
%         meta_obj_string_idx = find(covered_obj_string_idx_list == obj_string_idx); 
%         all_penalty_term_init = covered_obj_string_init_penalty{meta_obj_string_idx}; 
% 
%         % Compare penalty_term_init in resultFile to all_penalty_term_init
%             % Compare the 4th and 6th penalty 
%         resultFile_penalty_term_partial = [resultFile.penalty_term_init(4),resultFile.penalty_term_init(6)]; 
%         all_penalty_term_init_partial = [all_penalty_term_init(:,4),all_penalty_term_init(:,6)]; 
%         [~,row_idx] = ismember(resultFile_penalty_term_partial,all_penalty_term_init_partial,'rows');
% 
%         % Rename resultFile 
%         % new_resultFileName = sprintf('param_est_run_save/%s_iter_%d.mat',file_name(1:end-4),row_idx); 
%         % movefile(strcat('param_est_run_save/',file_name),new_resultFileName); 
% 
%         % end
%         % end
%     end
% end


%% Generate more objective string 
% target_file_name = 'qualitative_trend_obj_string.mat'; 
% target_file = load(target_file_name); 
% all_obj_string = target_file.all_obj_string; 
% original_obj_length = length(all_obj_string); 
%     % Going to add 
%     % crosstalk_ratio_transition_penalty,crosstalk_ratio_promotor_strength_penalty,log_baseline_data_penalty
% new_obj_string_list = {'log10(baseline_data_dev)','log10(crosstalk_ratio_dev)','log_baseline_data_penalty',...
%     'crosstalk_ratio_transition_penalty','crosstalk_ratio_promotor_strength_penalty','log10(other_penalty)','promotor_sensitivity_penalty',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + crosstalk_ratio_transition_penalty + 1000 * crosstalk_ratio_promotor_strength_penalty + 10000 * promotor_sensitivity_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + crosstalk_ratio_transition_penalty + 1000 * crosstalk_ratio_promotor_strength_penalty + 1000 * promotor_sensitivity_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + crosstalk_ratio_transition_penalty + 1000 * crosstalk_ratio_promotor_strength_penalty + 100 * promotor_sensitivity_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + 10 * crosstalk_ratio_transition_penalty + 1000 * crosstalk_ratio_promotor_strength_penalty + 10000 * promotor_sensitivity_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + 10 * crosstalk_ratio_transition_penalty + 1000 * crosstalk_ratio_promotor_strength_penalty + 1000 * promotor_sensitivity_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + 10 * crosstalk_ratio_transition_penalty + 1000 * crosstalk_ratio_promotor_strength_penalty + 100 * promotor_sensitivity_penalty + log10(other_penalty)'};
% updated_all_obj_string = [all_obj_string new_obj_string_list];
% all_obj_string = updated_all_obj_string; 
% save(target_file_name,'all_obj_string'); 

%% Test updated fit qualitative objective function flow 
    % Calculate the terms using an existing result file 

% mode = 'PE'; 
% num_conc = 7;
% num_promotor = 4; 
% 
%     % Load result file 
% resultFileName = 'param_est_run_save/20240510_param_est_run1234_objString_9.mat';
% resultFile = load(resultFileName); 
% problemObject = resultFile.problemObject;
% simFunction = resultFile.simFunction; 
% species_name_list = {problemObject.Model.Species.Name};
% experimental_grouped_data = problemObject.Data; 
% fitResult = resultFile.all_fitResults{resultFile.sort_idx(1),1}; 
% estimated_params = fitResult.estimated_params; 
% 
%     % Load objective string 
% obj_string_file_name = 'qualitative_trend_obj_string.mat';
% obj_string_file = load(obj_string_file_name);
% all_obj_string = obj_string_file.all_obj_string;
% obj_string = all_obj_string{9}; 
% 
%     % Define sensitivity option 
% rxn_oi_list = {'[CUTP:AGTP:RNAP:DNA pJ23105--utrGFP--sfGFP] -> [term_RNAP:DNA pJ23105--utrGFP--sfGFP] + [RNA utrGFP--sfGFP]'}; 
% param_oi_list = {'TXTL_PJ23105_RNAPbound_R'}; 
% modify_value_magnitude_list = -4:1:4; 
% modify_values = 10.^ modify_value_magnitude_list;
% modify_value_list = {modify_values}; 
% sensitivity_goal_list = {1}; 
% sensitivity_options = struct('rxn_oi_list',rxn_oi_list,'param_oi_list',param_oi_list,'modify_value_list',modify_value_list,...
%     'sensitivity_goal_list',sensitivity_goal_list);
% 
% [obj,baseline_data_dev,crosstalk_ratio_dev,log_baseline_data_penalty,crosstalk_ratio_transition_penalty,...
%     crosstalk_ratio_promotor_strength_penalty,promotor_sensitivity_penalty,other_penalty] = ...
%     wrapper_calculate_obj(estimated_params,obj_string,simFunction,problemObject,experimental_grouped_data,...
%     species_name_list,sensitivity_options,mode,num_conc,num_promotor);

    % Run an optimization
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%         'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%         'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%         'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'param_info/parameters_v2.xlsx';
% time = datetime; 
% save_path_stem = sprintf('param_est_run_save/2024%02d%02d_param_est_run%02d%02d',time.Month,time.Day,time.Hour,time.Minute); 
% obj_string_idx = 28;
% wrapper_fit_qualitative_characteristics(group_description,param_info_path,save_path_stem,obj_string_idx); 

%% Develop a workflow to check sensitivity of mRNA production to promotor binding affinity 
    % Adapted into calculate_sensitivity_penalty 
%     % mRNA production can be approximated by the 'TX elongation' step 
%     % Regardless of the fitted binding affinity (i.e. RNAP_bound_R), Set
%     % this value to be 1e-04:1e+04, solve the system, and record mRNA
%     % production rate? 
% 
%     % Take the TL-induced toxin on Ribosome mechanism fitting result as a
%     % case study 
% test_resultFileName = 'param_est_run_save/20240510_param_est_run1234_objString_9.mat';
% test_resultFile = load(test_resultFileName); 
% test_problemObject = test_resultFile.problemObject;
%     % Choose a fitResult to start with
% fitResult_oi = test_resultFile.all_fitResults{test_resultFile.sort_idx(1),1}; 
% estimated_params_name_list = {test_problemObject.Estimated.Name};
% estimated_params_value_list = fitResult_oi.estimated_params; 
% model = test_problemObject.Model; 
% species_name_list = {model.Species.Name};
% remove_idx = strcmp('protein sfGFP*',species_name_list); 
% additional_track_species_list = species_name_list(~remove_idx); 
% 
%     % Generate simData using simFunction 
% simFunction = create_simFun_from_problemObject(test_problemObject,additional_track_species_list); 
% dosing_information = create_dosing_info_from_problemObject(test_problemObject); 
% tStart = 0;
% tEnd = 21600; 
% [simulated_time,simulated_data] = simFunction(estimated_params_value_list,tEnd,dosing_information,tStart:tEnd);
% 
% plot_simulated(simulated_time,simulated_data,'baseline')
% plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% 
%     % Choose simData for sigma70 weak lowest plasmid concentration 
% simulated_time_single = simulated_time{85}; 
% simulated_data_single = simulated_data{85}; 
% simulated_data_name_list = [{'protein sfGFP*'},additional_track_species_list];
% 
% [rxn_name_list,rxn_time_course] = get_reaction_time_course(model,estimated_params_name_list,estimated_params_value_list,...
%     simulated_time_single,simulated_data_single,simulated_data_name_list);
% 
%     % Get the control mRNA production rate 
%         % Note: need to validate this 
% rxn_oi = '[CUTP:AGTP:RNAP:DNA pJ23105--utrGFP--sfGFP] -> [term_RNAP:DNA pJ23105--utrGFP--sfGFP] + [RNA utrGFP--sfGFP]';
% rxn_idx = strcmp(rxn_name_list,rxn_oi);
% rxn_oi_time_course_ctrl = rxn_time_course(:,rxn_idx); 
% % legend_names = {'Control'};
% % figure;
% % plot(simulated_time_single,rxn_oi_time_course_ctrl,'LineWidth',1.5,'Color','k')
% % hold on 
% 
%     % Vary RNAP_bound_R and record the mRNA production rate 
%         % Find index of RNAP_bound_J23105_R
% modify_param_oi = 'TXTL_PJ23105_RNAPbound_R';
% modify_param_idx = strcmp(estimated_params_name_list,modify_param_oi); 
% 
%     % This can be recorded as fold changes of max RNA production rate across these binding affinities 
% modify_value_magnitude_list = -4:1:4; 
% modify_value_list = 10.^ modify_value_magnitude_list;
% all_max_RNA_prod = nan(length(modify_value_list),1); % Record max RNA production rate 
% 
% for modify_idx = 1:length(modify_value_list)
%     modify_value = modify_value_list(modify_idx); 
%     modified_estimated_params_value_list = estimated_params_value_list; 
%     modified_estimated_params_value_list(modify_param_idx) = modify_value; 
% 
%     [simulated_time,simulated_data] = simFunction(modified_estimated_params_value_list,tEnd,dosing_information,tStart:tEnd);
% 
%     simulated_time_single = simulated_time{85}; 
%     simulated_data_single = simulated_data{85}; 
%     simulated_data_name_list = [{'protein sfGFP*'},additional_track_species_list];
% 
%     [rxn_name_list,rxn_time_course] = get_reaction_time_course(model,estimated_params_name_list,modified_estimated_params_value_list,...
%     simulated_time_single,simulated_data_single,simulated_data_name_list);
%     rxn_oi_time_course = rxn_time_course(:,rxn_idx); 
% 
%     all_max_RNA_prod(modify_idx) = max(rxn_oi_time_course); 
% 
%     % plot(simulated_time_single,rxn_oi_time_course,'LineWidth',1.5)
%     % legend_names{end+1} = num2str(modify_value); 
% end
% max_fold_change_mRNA_prod = all_max_RNA_prod(1) / all_max_RNA_prod(end); 
% % legend(legend_names)
%% Function to calculate time-course reaction rates from concentration profiles 
% 
% model = problemObject.Model;
% simData = fitted_data(1); 
% 
% 
% [rxn_name_list,rxn_time_course] = get_reaction_time_course(model,fitResult,simData); 
% % Get model reactions, parameters, and fitted parameters 
% model_reaction_obj_list = problemObject.Model.Reactions; 
% model_parameter_obj_list = problemObject.Model.Parameters; 
% estimated_parameter_name_list = fitResult_oi.ParameterEstimates.Name; 
% estimated_parameter_value_list = [fitResult_oi.ParameterEstimates.Estimate];
% 
% % Get tracked species name and time-course
% test_fitted_data = fitted_data(1,1).Data; 
% test_fitted_data_name_list = fitted_data(1,1).DataNames; 
% test_fitted_data_time_list = fitted_data(1,1).Time; 
% 
% % Assign time-course to each species 
% for species_idx = 1:length(test_fitted_data_name_list)
%     species_name = test_fitted_data_name_list{species_idx};
%         % Replace some characters to avoid assigning error
%     if strcmp(species_name,'protein sfGFP*')
%         species_name = 'protein_sfGFP_m';
%     end
%     species_name = strrep(species_name,' ','');
%     species_name = strrep(species_name,'--',''); 
%     species_name = strrep(species_name,':',''); 
% 
%     eval(sprintf('%s = test_fitted_data(:,species_idx);',species_name)); 
% 
% end
% 
% % Assign parameter values 
%     % Existing parameters in problemObject
% for param_idx_1 = 1:length(model_parameter_obj_list)
%     parameter_obj = model_parameter_obj_list(param_idx_1); 
%     parameter_name = parameter_obj.Name; 
%     eval(sprintf('%s = parameter_obj.Value;',parameter_name)); 
% end
%     % Update these based on fitResult
% for param_idx_2 = 1:length(estimated_parameter_name_list)
%     parameter_name = estimated_parameter_name_list{param_idx_2}; 
%     eval(sprintf('%s = estimated_parameter_value_list(param_idx_2);',parameter_name)); 
% end
% 
% % Preassign a matrix for reaction rate calculation 
% rxn_time_course = nan(length(test_fitted_data_time_list),length(model_reaction_obj_list)); 
% 
% for rxn_idx = 1:length(model_reaction_obj_list)
%     reaction_obj = model_reaction_obj_list(rxn_idx); 
%     rxn_rate_string = reaction_obj.ReactionRate; 
%         % Replace strings accordingly compared to previous assignments 
%     rxn_rate_string = strrep(rxn_rate_string,' ','');
%     rxn_rate_string = strrep(rxn_rate_string,'--','');
%     rxn_rate_string = strrep(rxn_rate_string,':','');
%     rxn_rate_string = strrep(rxn_rate_string,'protein_sfGFP*','protein_sfGFP_m');
%     rxn_rate_string = strrep(rxn_rate_string,'[','');
%     rxn_rate_string = strrep(rxn_rate_string,']','');
%         % Get .* for correct dimension
%     rxn_rate_string = strrep(rxn_rate_string,'*','.*'); 
% 
%     rxn_time_course(:,rxn_idx) = eval(rxn_rate_string); 
% 
% 
% end


%% Check different toxin mechanism 
% 
%     % Result file 
% ctrl_TX_gen_on_RNAP_resultFileName ='param_est_run_save/20231101_param_est_run0929_RNAP_toxin_mech.mat';
% TX_gen_on_Ribo_resultFileName = 'param_est_run_save/20230125_param_est_run1111_TXgen_Ribo.mat';
% TL_gen_on_RNAP_resultFileName = 'param_est_run_save/20230125_param_est_run2230_TLgen_RNAP.mat';
% TL_gen_on_Ribo_resultFileName = 'param_est_run_save/20230125_param_est_run1326_TLgen_Ribo.mat';
% 
% % % Plot out T7 baseline & see whether decrease in expression is captured 
% %     % Define result analysis function inputs 
% % metric_name = 'LogLikelihood'; 
% % plot_fitting = true;
% % plot_param_dist = false;
% % 
% %     % Control (TX generated toxin poisons RNAP) 
% % ctrl_TX_gen_on_RNAP_resultFile = load(ctrl_TX_gen_on_RNAP_resultFileName);
% % % wrapper_analyze_fit_result_v2(ctrl_TX_gen_on_RNAP_resultFile.all_fitResults,ctrl_TX_gen_on_RNAP_resultFile.problemObject,...
% % %     ctrl_TX_gen_on_RNAP_resultFileName,metric_name,plot_fitting,plot_param_dist); 
% % ctrl_TX_gen_on_RNAP_metric_summary = ctrl_TX_gen_on_RNAP_resultFile.metric_summary; 
% % [~,min_idx] = min(ctrl_TX_gen_on_RNAP_metric_summary); 
% % ctrl_TX_gen_on_RNAP_fitResult_oi = ctrl_TX_gen_on_RNAP_resultFile.all_fitResults{min_idx,1}; 
% % [Time,Data]
% % plot_simulated(Time,Data,option,exp_Time,exp_Data); 
% % 
% %     % TX generated toxin poisons Ribosome 
% % TX_gen_on_Ribo_resultFile = load(TX_gen_on_Ribo_resultFileName); 
% % % wrapper_analyze_fit_result_v2(TX_gen_on_Ribo_resultFile.all_fitResults,TX_gen_on_Ribo_resultFile.problemObject,...
% % %     TX_gen_on_Ribo_resultFileName,metric_name,plot_fitting,plot_param_dist); 
% % 
% %     % TL generated toxin poisons RNAP
% % TL_gen_on_RNAP_resultFile = load(TL_gen_on_RNAP_resultFileName);
% % % wrapper_analyze_fit_result_v2(TL_gen_on_RNAP_resultFile.all_fitResults,TL_gen_on_RNAP_resultFile.problemObject,...
% % %     TL_gen_on_RNAP_resultFileName,metric_name,plot_fitting,plot_param_dist); 
% % 
% %     % TL generated toxin poisons Ribosome 
% % TL_gen_on_Ribo_resultFile = load(TL_gen_on_Ribo_resultFileName);
% % % wrapper_analyze_fit_result_v2(TL_gen_on_Ribo_resultFile.all_fitResults,TL_gen_on_Ribo_resultFile.problemObject,...
% % %     TL_gen_on_Ribo_resultFileName,metric_name,plot_fitting,plot_param_dist); 
% 
% load(TL_gen_on_RNAP_resultFileName); 
% 
% figure;
% for iter = 1:length(all_fitResults)
%     fitResult_oi = all_fitResults{iter}; 
%     fitted_data = fitted(fitResult_oi); 
% 
%     [Time,Data] = process_crosstalk_data_from_source(fitted_data,'simulated_data');
%     % [Time,Data] = process_crosstalk_data_from_source(problemObject.Data,'grouped_data');
%     subplot(7,7,iter)
%     for conc_idx = 1:length(Data)
%         plot(Time{conc_idx}./3600,Data{conc_idx},'LineWidth',1.5)
%         hold on 
%     end
%     sgtitle('TL-gen toxin on RNAP')
%     % xlabel('Time (hr)','FontSize',14)
%     % ylabel('Concentration (nM)','FontSize',14)
%     % legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')
%     set(gca,'FontSize',14)
% end


%% Add new objective strings
% target_file_name = 'qualitative_trend_obj_string.mat'; 
% target_file = load(target_file_name); 
% all_obj_string = target_file.all_obj_string; 
% original_obj_length = length(all_obj_string); 
%     % Going to add 
%     % crosstalk_ratio_transition_penalty,crosstalk_ratio_promotor_strength_penalty,log_baseline_data_penalty
% new_obj_string_list = {'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + crosstalk_ratio_transition_penalty + crosstalk_ratio_promotor_strength_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + crosstalk_ratio_transition_penalty + 10^crosstalk_ratio_promotor_strength_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + crosstalk_ratio_transition_penalty + 10^(-1 * crosstalk_ratio_promotor_strength_penalty) + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + crosstalk_ratio_transition_penalty + 1000 * crosstalk_ratio_promotor_strength_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + 10 * crosstalk_ratio_transition_penalty + 10^crosstalk_ratio_promotor_strength_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + 10 * crosstalk_ratio_transition_penalty + 10^(-1 * crosstalk_ratio_promotor_strength_penalty) + log10(other_penalty)',...
%     'log10(baseline_data_dev) + log10(crosstalk_ratio_dev) + log_baseline_data_penalty + 10 * crosstalk_ratio_transition_penalty + 1000 * crosstalk_ratio_promotor_strength_penalty + log10(other_penalty)'};
% updated_all_obj_string = [all_obj_string new_obj_string_list];
% all_obj_string = updated_all_obj_string; 
% save(target_file_name,'all_obj_string'); 

%% Take fitted results and calculate the penalty terms 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'benchmark_model/parameters_test_PE.xlsx';
% problemObject = setProblemObject_PE(group_description,param_info_path,[]);
% tStart = 0; 
% tEnd = 21600; % Modify this later to draw info from data  
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% all_obj_string_file = load('qualitative_trend_obj_string.mat'); 
% all_obj_string = all_obj_string_file.all_obj_string; 
% 
% obj_string = all_obj_string{8}; 
% result_fileName = 'param_est_run_save/20240125_param_est_run1110_objString_8.mat';
% result_file = load(result_fileName); 
% simFunction = result_file.simFunction; 
% experimental_data = problemObject.Data;
% species_name_list = {problemObject.Model.Species.Name};
% mode = 'PE'; 
% num_conc = 7; 
% num_promotor = 4; 
% opt_fitResult = result_file.all_fitResults{result_file.sort_idx(2)}; 
% fitted_params = opt_fitResult.estimated_params; 
% 
% obj_val = wrapper_calculate_obj(fitted_params,obj_string,simFunction,problemObject,experimental_data,...
%     species_name_list,mode,num_conc,num_promotor); 
%% Check results from fitting various objective weights 
 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'benchmark_model/parameters_test_PE.xlsx';
% problemObject = setProblemObject_PE(group_description,param_info_path,[]);
% tStart = 0; 
% tEnd = 21600; % Modify this later to draw info from data  
% % 
% % all_obj_string_file = load('qualitative_trend_obj_string.mat'); 
% % all_obj_string = all_obj_string_file.all_obj_string; 
% 
% time_stamps = {'1108','1109','1110'};
% time_stamps = {'1234','1246'}; 
% for obj_string_idx = 9:2:13
%     if obj_string_idx >= 1 
%         obj_string = all_obj_string{obj_string_idx}; 
%     else
%         obj_string = 'log10(baseline_data_dev) + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + crosstalk_ratio_penalty + log10(other_penalty)'; 
%     end
%     for time_stamp_idx = 1:length(time_stamps)
%         time_stamp = time_stamps{time_stamp_idx}; 
%         % resultFile_objString = sprintf('param_est_run_save/20240125_param_est_run%s_objString_%d.mat',time_stamp,obj_string_idx); 
%         resultFile_objString = sprintf('param_est_run_save/20240510_param_est_run%s_objString_%d.mat',time_stamp,obj_string_idx); 
%         if exist(resultFile_objString,'file')
%             resultFile = load(resultFile_objString); 
%             break
%         end
%     end
% 
%     % Check baseline and crosstalk ratio 
%     problemObject = resultFile.problemObject;
%     dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
%     simFunction = resultFile.simFunction; 
%     all_fitResults = resultFile.all_fitResults; 
%     all_fval = nan(length(all_fitResults),1); 
%     for iter = 1:length(all_fitResults)
%         fitResult = all_fitResults{iter,1}; 
%         all_fval(iter,1) = fitResult.fval; 
%     end
% 
%     [sorted_fval,sort_idx] = sort(all_fval,'ascend'); 
% 
%     for idx = 1:5
% 
%         fitResult = all_fitResults{sort_idx(idx),1}; 
%         estimated_params = fitResult.estimated_params; 
%         [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
%         plot_simulated(simulated_time,simulated_data,'baseline')
%         plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
%         sgtitle(sprintf('%s \n fval = %.2f',obj_string,sorted_fval(idx)))
% 
%     end
%     save(resultFile_objString,'all_fval','sorted_fval','sort_idx','-append')
% 
%     % Track penalty value change 
%         % Break obj string down into segments for eval statements 
%     separate_idx_list = strfind(obj_string,'+'); 
%     penalty_term_str_list = cell(length(resultFile.penalty_term_label),1);
%     for sep_idx = 1:length(separate_idx_list) + 1 
%         if isequal(sep_idx,1)
%             start_idx = 1; 
%         else
%             start_idx = end_idx; 
%         end
%         if isequal(sep_idx,length(separate_idx_list) + 1)
%             end_idx = length(obj_string); 
%         else
%             end_idx = separate_idx_list(sep_idx); 
%         end
%         term_oi = obj_string(start_idx:end_idx);
%         if contains(term_oi,'^')
%             term_oi = strrep(term_oi,'^','.^');
%         end
%         penalty_term_str_list{sep_idx} = strrep(term_oi,'+','');
%     end
% 
%         % Use a stacked bar plot in each case 
%     figure; 
%     for iter = 1:size(resultFile.all_penalty_term_final,1)
%         penalty_term_init = resultFile.all_penalty_term_init(iter,:);
%         penalty_term_final = resultFile.all_penalty_term_final(iter,:); 
%         penalties = [penalty_term_init;penalty_term_final];
%         modified_penalties = penalties; 
%         for penalty_idx = 1:length(penalty_term_str_list)
%             penalty_term_label = resultFile.penalty_term_label{penalty_idx}; 
%             eval(sprintf('%s = penalties(:,penalty_idx);',penalty_term_label)); 
%             eval(sprintf('modified_penalty = %s;',penalty_term_str_list{penalty_idx})); 
%             modified_penalties(:,penalty_idx) = modified_penalty; 
%         end
%         subplot(7,7,iter)
%         bar(modified_penalties,'stacked')
%     end
%     legend(strrep(resultFile.penalty_term_label,'_',' '))
%     sgtitle(sprintf('objective string %d',obj_string_idx))
% 
% end


%% Temp plotting for group meeting 
% crosstalk_ratios = calculate_crosstalk_ratio_v2(Time,Data,length(conc_vec),4,'PE'); 
% figure; 
% for prom_idx = 1:4
%     subplot(2,2,prom_idx)
%     promotor_name = promotor_name_list{prom_idx}; 
%     prom_crosstalk_ratio = crosstalk_ratios(:,prom_idx); 
% 
%     if ~isequal(prom_idx,4)
%         scatter(conc_vec(3:end),prom_crosstalk_ratio{1}(3:end),'MarkerFaceColor',all_colors{2})
%         hold on 
%         scatter(conc_vec(3:end),prom_crosstalk_ratio{2}(3:end),'MarkerFaceColor',all_colors{4})
%         scatter(conc_vec(3:end),prom_crosstalk_ratio{3}(3:end),'MarkerFaceColor',all_colors{3})
%         plot(conc_vec(3:end),prom_crosstalk_ratio{1}(3:end),'Color',all_colors{2},'LineWidth',1.5)
%         hold on 
%         plot(conc_vec(3:end),prom_crosstalk_ratio{2}(3:end),'Color',all_colors{4},'LineWidth',1.5)
%         plot(conc_vec(3:end),prom_crosstalk_ratio{3}(3:end),'Color',all_colors{3},'LineWidth',1.5)
%     else
%         scatter(conc_vec(5:end),prom_crosstalk_ratio{1}(5:end),'MarkerFaceColor',all_colors{2})
%         hold on 
%         scatter(conc_vec(5:end),prom_crosstalk_ratio{2}(5:end),'MarkerFaceColor',all_colors{4})
%         scatter(conc_vec(5:end),prom_crosstalk_ratio{3}(5:end),'MarkerFaceColor',all_colors{3})
%         plot(conc_vec(5:end),prom_crosstalk_ratio{1}(5:end),'Color',all_colors{2},'LineWidth',1.5)
%         hold on 
%         plot(conc_vec(5:end),prom_crosstalk_ratio{2}(5:end),'Color',all_colors{4},'LineWidth',1.5)
%         plot(conc_vec(5:end),prom_crosstalk_ratio{3}(5:end),'Color',all_colors{3},'LineWidth',1.5)
%     end
%     title(strrep(promotor_name,'_',' '))
%     legend('empty','empty T7','empty sigma70')
% end

%% Generate objective strings for qualitative trends 
% Before we generate, let's check the results from large-scale fitting 
% large_scale_fitting_resultFile = load('param_est_run_save/20231106_large_scale_run.mat'); 
% 
% all_obj_string = {'log10(baseline_data_dev) + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + 10 * crosstalk_ratio_penalty + log10(other_penalty)',...
%     'log10(baseline_data_dev) + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + 100 * crosstalk_ratio_penalty + log10(other_penalty)',...
%     'baseline_data_dev/10^11 + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + crosstalk_ratio_penalty + log10(other_penalty)',...
%     'baseline_data_dev/10^11 + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + 10 * crosstalk_ratio_penalty + log10(other_penalty)',...
%     'baseline_data_dev/10^11 + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + 100 * crosstalk_ratio_penalty + log10(other_penalty)',...
%     'baseline_data_dev/10^12 + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + crosstalk_ratio_penalty + log10(other_penalty)',...
%     'baseline_data_dev/10^12 + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + 10 * crosstalk_ratio_penalty + log10(other_penalty)',...
%     'baseline_data_dev/10^12 + crosstalk_ratio_dev/10 + log10(baseline_data_penalty) + 100 * crosstalk_ratio_penalty + log10(other_penalty)'}; 
% save('qualitative_trend_obj_string.mat','all_obj_string'); 
%% Check results for qualitative trend fitting 
% resultFile = load('param_est_run_save/20230116_param_est_run1020.mat'); 
% simFunction = resultFile.simFunction; 
% 
% 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%         'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%         'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%         'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'benchmark_model/parameters_test_PE.xlsx';
% problemObject = setProblemObject_PE(group_description,param_info_path,[]);
% tStart = 0; 
% tEnd = 21600; % Modify this later to draw info from data  
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
% for idx = 1:length(resultFile.all_fitResults)
%     fitResult = resultFile.all_fitResults{idx,1}; 
%     estimated_params = fitResult.estimated_params; 
%     [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
%     plot_simulated(simulated_time,simulated_data,'baseline')
% 
% end


%% Choose a run as a 'toy model' to test for large positive crosstalk ratio in weak promotors 
% %     % Thinking - pick a set where crosstalk ratio is sensitive to promotor
% %     % binding affinities? This could be hard to implement though - maybe
% %     % just choosing a set with diverse crosstalk ratios for promotors would
% %     % be okay? (e.g. big difference between T7 strong and weak? ) 
% % % load param_est_run_save/20231106_large_scale_run.mat
% % % data_table = updated_optimization_run_table; 
% % % T7_diff = nan(height(data_table),7); 
% % % sigma70_diff = nan(height(data_table),7); 
% % % pos_flag_array = nan(height(data_table),1); 
% % % 
% % %     % Choose a fitted result as a starting point to work with 
% % % for row_idx = 1:height(data_table)
% % %     crosstalk_ratios = data_table{row_idx,'CrosstalkRatios'}{1,1};
% % %     [T7_strong_crosstalk_ratio,T7_weak_crosstalk_ratio,sigma70_strong_crosstalk_ratio,sigma70_weak_crosstalk_ratio] = ...
% % %         unpack_crosstalk_ratios(crosstalk_ratios,promotor_name_list); 
% % %     T7_diff(row_idx,:) = T7_strong_crosstalk_ratio(:,1)' - T7_weak_crosstalk_ratio(:,1)'; 
% % %     sigma70_diff(row_idx,:) = sigma70_strong_crosstalk_ratio(:,1)' - sigma70_weak_crosstalk_ratio(:,1)'; 
% % %     pos_flag_array(row_idx) = any(sigma70_weak_crosstalk_ratio(1:3,1) > 1.1); 
% % % end
% % % % [sorted_sigma70_diff,sort_idx] = sort(abs(sigma70_diff(:,1)),'descend'); 
% % % % row_oi = data_table(sort_idx(1),:);
% % % eligible_row_list = find(pos_flag_array); 
% % % pos_crosstalk_sigma70_diff = sigma70_diff(eligible_row_list,:); 
% % % row_oi = data_table(eligible_row_list(12),:);
% % 
% %     % Load the fitting result and get simFunction 
% % resultFile = load('param_est_run_save/20231102_param_est_run1534_96.mat'); 
% fitResult_oi = resultFile.all_fitResults{14,1}; 
% probObject = resultFile.problemObject; 
% 
% 
%     % Additional species to track
% species_name_list = {probObject.Model.Species.Name};
% mRNA_idx_list = contains(species_name_list,'RNA ');
% mRNA_species = species_name_list(mRNA_idx_list); 
% DNA_idx_list = contains(species_name_list,'DNA '); 
% DNA_species = species_name_list(DNA_idx_list); 
% additional_track_species = [{'RNAP','RNase','Ribo','t7RNAP'},mRNA_species,DNA_species]; 
% simFunction = create_simFun_from_problemObject(probObject,additional_track_species); 
% % simFunction = create_simFun_from_problemObject(probObject); 
% dosing_information = create_dosing_info_from_problemObject(probObject); 
% param_names = {probObject.Estimated.Name};
% ori_param_val = [fitResult_oi.ParameterEstimates.Estimate];
% 
%     % Modify parameters 
% tStart = 0;
% tEnd = 21600; 
% % modify_param_name_list = {'TXTL_PkanR_RNAPbound_R','TXTL_PJ23105_RNAPbound_R','tx_capacity_param','TL_elong_glob',...
% %     'TX_elong_glob','TX_elong_glob_T7','TXTL_RNAdeg_R','RNase_0',...
% %     'TXTL_RNAdeg_R_sfGFP','TXTL_RNAdeg_R','TXTL_RNAdeg_kc_sfGFP','TXTL_RNAdeg_kc_kanR'}; 
% % modify_param_val_list = [2000000,20000000,20,5,...
% %     60,40,0.000001,3,...
% %     0.001,0.000001,0.000001,0.000001];
% modify_param_name_list = {'TXTL_PkanR_RNAPbound_R','TXTL_PJ23105_RNAPbound_R','tx_capacity_param','TL_elong_glob',...
%     'TX_elong_glob','TX_elong_glob_T7','TXTL_RNAdeg_R','RNase_0',...
%     'TXTL_RNAdeg_R_sfGFP','TXTL_RNAdeg_R','TXTL_RNAdeg_kc_sfGFP','TXTL_RNAdeg_kc_kanR',...
%     'AGTPdeg_time'}; 
% modify_param_val_list = [2000000,20000000,20,5,...
%     60,40,0.000001,30,...
%     0.001,0.000001,0.0001,0.0001,...
%     1800];
% % modify_param_name_list = {'TXTL_PkanR_RNAPbound_R','TXTL_PJ23105_RNAPbound_R','tx_capacity_param','TL_elong_glob',...
% %     'TX_elong_glob','TX_elong_glob_T7','RNase_0',...
% %     'TXTL_RNAdeg_R_sfGFP','TXTL_RNAdeg_R','TXTL_RNAdeg_kc_sfGFP','TXTL_RNAdeg_kc_kanR','TXTL_RNAdeg_kc',...
% %     'AGTPdeg_time','RNAP_0'}; 
% % modify_param_val_list = [2000000,1000000,20,5,...
% %     60,40,15,...
% %     10,0.01,0.001,0.001,1,...
% %     1800,800];
% modified_params = modify_params(param_names,ori_param_val,modify_param_name_list,modify_param_val_list); 
% [Time,Data] = simFunction(modified_params',tEnd,dosing_information,tStart:tEnd);
% 
%     % % One more sanity check to make sure results from fitted & simFunction
%     % % are the same 
%     % simulated_data = fitted(fitResult_oi); 
%     % [simTime,simData] = process_crosstalk_data_from_source(simulated_data,'simulated_data'); 
%     % 
%     % 
%     % figure; 
%     % for idx = 1:81% length(Data)
%     %     simFunc_Time = Time{idx,1}; 
%     %     simFunc_Data = Data{idx,1}; 
%     %     fittedFunc_Time = simTime{idx,1}; 
%     %     fittedFunc_Data = simData{idx,1}; 
%     %     subplot(9,9,idx)
%     % 
%     %     plot(fittedFunc_Time,fittedFunc_Data,'r','LineWidth',1.5)
%     %     hold on 
%     %     plot(simFunc_Time,simFunc_Data,'k','LineWidth',1.5);
%     % end
% 
%     % Track the impact of changing promotor binding affinity on expression
%     % & positive crosstalk ratio 
% 
% % Plot out (1) baseline expression T7 strong, T7 weak, sigma70 strong,
% % sigma70 weak (2) Crosstalk ratio
% 
% plot_simulated(Time,Data,'crosstalk_ratio')
% plot_simulated(Time,Data,'baseline')
% 
% % At this point, probably a good idea to check the produced mRNA species
% % for each case - this is to see how much reporter vs. distraction are
% % produced and how much are bound 
% T7_strong_data_idx_list = [1,8,15,22];
% sigma70_strong_data_idx_list = [57,64,71,78];
% sigma70_weak_data_idx_list = [85,92,99,106];
% 
% % % Check other reporter plasmid concentrations 
% % T7_strong_data_idx_list = T7_strong_data_idx_list + 3; 
% % sigma70_strong_data_idx_list = sigma70_strong_data_idx_list + 3; 
% % sigma70_weak_data_idx_list = sigma70_weak_data_idx_list + 3; 
% 
% Data_species_list = [{'sfGFP*'},additional_track_species];
% T7_strong_Time = Time(T7_strong_data_idx_list); 
% T7_strong_Data = Data(T7_strong_data_idx_list); 
% figure_title = 'T7 strong 0.5nM';
% 
% plot_mRNA(T7_strong_Time,T7_strong_Data,Data_species_list,{'RNase'},figure_title)
% plot_mRNA(Time(sigma70_strong_data_idx_list),Data(sigma70_strong_data_idx_list),Data_species_list,{'RNase'},'sigma70 strong 0.5nM')
% plot_mRNA(Time(sigma70_weak_data_idx_list),Data(sigma70_weak_data_idx_list),Data_species_list,{'RNase'},'sigma70 weak 0.5nM')
% 
% plot_DNA(T7_strong_Time,T7_strong_Data,Data_species_list,{'t7RNAP','RNAP'},figure_title)
% plot_DNA(Time(sigma70_strong_data_idx_list),Data(sigma70_strong_data_idx_list),Data_species_list,{'t7RNAP','RNAP'},'sigma70 strong 0.5nM')
% plot_DNA(Time(sigma70_weak_data_idx_list),Data(sigma70_weak_data_idx_list),Data_species_list,{'t7RNAP','RNAP'},'sigma70 weak 0.5nM')
%% Add missing large-scale fitting result to structure 

%     % Process raw rerun data 
% result_fileName_stem = 'param_est_run_save/20231108_param_est_run0936'; 
% rerun_idx_list = [5,67,86];
% for idx = 1:length(rerun_idx_list)
%     rerun_idx = rerun_idx_list(idx); 
%     result_fileName_suffix = sprintf('_%d.mat',rerun_idx); 
%     result_fileName = strcat(result_fileName_stem,result_fileName_suffix); 
%     result_file = load(result_fileName); 
% 
%     all_fitResults = result_file.all_fitResults; 
%     problemObject = result_file.problemObject; 
%     fit_result_path = result_fileName(1:end-4); 
%     metric_name = 'SSE'; 
%     plot_fitting = false;
%     plot_param_dist = false; 
%     wrapper_analyze_fit_result_v2(all_fitResults,problemObject,fit_result_path,metric_name,plot_fitting,plot_param_dist) 
% 
% end

%     % Add the rerun data to stored result structure 
% load param_est_run_save/20231106_large_scale_run.mat
% updated_optimization_run_table = optimization_run_table; 
% rerun_idx_list = [5,67,86];
% num_conc = 7;
% num_promotor = 4; 
% for idx = 1:length(rerun_idx_list)
%     rerun_idx = rerun_idx_list(idx); 
%     result_fileName = sprintf('param_est_run_save/20231108_param_est_run0936_%d.mat',rerun_idx); 
%     result_file = load(result_fileName); 
%     for iter = 1:length(result_file.all_fitResults)
%         fitResult = result_file.all_fitResults{iter};
%         all_simData = fitted(fitResult); 
%         [Time,Data] = process_crosstalk_data_from_source(all_simData,'simulated_data'); 
%         crosstalk_ratio = calculate_crosstalk_ratio_v2(Time,Data,num_conc,num_promotor,'PE'); 
%         label = sprintf('run%d_iter%d',rerun_idx,iter); 
%         SSE = fitResult.SSE; 
%         LogLikelihood = fitResult.LogLikelihood; 
%         ParameterEstimates = [fitResult.ParameterEstimates.Estimate]; 
%         updated_optimization_run_table{(rerun_idx - 1) * 48 + iter,'SSE'} = SSE; 
%         updated_optimization_run_table{(rerun_idx - 1) * 48 + iter,'LogLikelihood'} = LogLikelihood; 
%         updated_optimization_run_table{(rerun_idx - 1) * 48 + iter,'CrosstalkRatios'} = {crosstalk_ratio}; 
%         updated_optimization_run_table{(rerun_idx - 1) * 48 + iter,'ParameterEstimates'} = {ParameterEstimates}; 
%     end
% end
% 
% save('param_est_run_save/20231106_large_scale_run.mat','updated_optimization_run_table'); 

% load('param_est_run_save/20231005_large_scale_run_1002.mat'); 
%     % # Load a result file 
% test_file = load('param_est_run_save/20231002_param_est_run1009_1.mat'); 
% param_name_list = {test_file.problemObject.Estimated.Name};
% 
% all_SSE = [optimization_run_table.SSE]; 
% all_estimated_parameters = [optimization_run_table.ParameterEstimates{:}]';
% 
% [sorted_all_SSE,sort_idx] = sort(all_SSE,'ascend','MissingPlacement','last'); 
% sorted_all_estimated_parameters = all_estimated_parameters(sort_idx,:); 
% 
% top10_perc_estimated_parameters = sorted_all_estimated_parameters(1:480,:); 
% [~,top10_perc_estimated_parameters_rank] = sort(top10_perc_estimated_parameters);
% 
% param_corr_mat = corrcoef(top10_perc_estimated_parameters_rank); 
% distance = 1 - abs(param_corr_mat); 
% Z = linkage(distance,'ward'); 
% cluster_labels = cluster(Z,'cutoff',1.5,'Criterion','distance'); 
% 
% figure;
% param_name_list_disp = strrep(param_name_list,'_',' '); 
% dendrogram(Z,'Labels',param_name_list_disp)
% 
% cluster_num = unique(cluster_labels);
% clustered_params = []; 
% clustered_param_name_list = {}; 
% num_in_cluster = nan(length(cluster_num),1); 
% for cluster_idx = 1:length(cluster_num)
%     param_idx_in_cluster = find(cluster_labels == cluster_idx); 
%     num_in_cluster(cluster_idx,1) = length(find(param_idx_in_cluster));
%     clustered_params = [clustered_params top10_perc_estimated_parameters_rank(:,param_idx_in_cluster)];
%     clustered_param_name_list = [clustered_param_name_list param_name_list(param_idx_in_cluster)]; 
% 
% end
% clustered_param_corr_mat = corrcoef(clustered_params); 
% 
% for param_idx = 1:size(clustered_param_corr_mat,1)
%     clustered_param_corr_mat(param_idx,param_idx) = nan;
% end
% figure; 
% clustered_param_name_list_disp = strrep(clustered_param_name_list,'_',' '); 
% h = heatmap(clustered_param_name_list_disp',clustered_param_name_list_disp',abs(clustered_param_corr_mat),'MissingDataColor','w');


%% Test parameter fitting with and without kanR TX/TL 
% resultFile_wo_kanR = load('param_est_run_save/20231129_param_est_run1325_mRNA_prod_deg.mat');
% resultFile_w_kanR_TX = load('param_est_run_save/20231127_param_est_run1707_mRNA_prod_deg_TXkanR.mat');
% resultFile_w_kanR_TXTL = load('param_est_run_save/20231127_param_est_run1728_mRNA_prod_deg_TXTLkanR.mat');
% 
% fitResultPath_wo_kanR = 'param_est_run_save/20231129_param_est_run1325_mRNA_prod_deg_analysis.mat';
% fitResultPath_w_kanR_TX = 'param_est_run_save/20231127_param_est_run1707_mRNA_prod_deg_TXkanR_analysis.mat';
% fitResultPath_w_kanR_TXTL = 'param_est_run_save/20231127_param_est_run1728_mRNA_prod_deg_TXTLkanR_analysis.mat';
% 
% metric_name = 'SSE';
% plot_fitting = false; 
% plot_param_dist = false; 
% 
% wrapper_analyze_fit_result_v2(resultFile_wo_kanR.all_fitResults,resultFile_wo_kanR.problemObject,...
%     fitResultPath_wo_kanR,metric_name,plot_fitting,plot_param_dist) 
% wrapper_analyze_fit_result_v2(resultFile_w_kanR_TX.all_fitResults,resultFile_w_kanR_TX.problemObject,...
%     fitResultPath_w_kanR_TX,metric_name,plot_fitting,plot_param_dist) 
% wrapper_analyze_fit_result_v2(resultFile_w_kanR_TXTL.all_fitResults,resultFile_w_kanR_TXTL.problemObject,...
%     fitResultPath_w_kanR_TXTL,metric_name,plot_fitting,plot_param_dist) 

% resultFile_wo_kanR = load('param_est_run_save/20231129_param_est_run1325_mRNA_prod_deg_analysis.mat');
% resultFile_w_kanR_TX = load('param_est_run_save/20231127_param_est_run1707_mRNA_prod_deg_TXkanR_analysis.mat');
% resultFile_w_kanR_TXTL = load('param_est_run_save/20231127_param_est_run1728_mRNA_prod_deg_TXTLkanR_analysis.mat');

%     % Compare the error metrics 
%           % These are all fitting to the same piece of data and the error
%           % metrics should be directly comparable 
% % Compare distribution of SSE 
% figure; 
% subplot(1,3,1)
% histogram(log10(resultFile_wo_kanR.metric_summary))
% xlabel('log(SSE)')
% ylabel('Occurance')
% title('without kanR')
% subplot(1,3,2)
% histogram(log10(resultFile_w_kanR_TX.metric_summary))
% title('with kanR TX')
% xlabel('log(SSE)')
% subplot(1,3,3)
% histogram(log10(resultFile_w_kanR_TXTL.metric_summary))
% title('with kanR TXTL')
% xlabel('log(SSE)')
% sgtitle('Distribution of fitting SSE')
% 
% % Compare min SSE 
% figure; 
% xlabels = categorical({'without kanR','with kanR TX','with kanR TXTL'}); 
% xlabels = reordercats(xlabels,{'without kanR','with kanR TX','with kanR TXTL'}); 
% bar(xlabels,[min(resultFile_wo_kanR.metric_summary),min(resultFile_w_kanR_TX.metric_summary),min(resultFile_w_kanR_TXTL.metric_summary)])
% title('min SSE')
% 
%     % Compare distribution of common estimated parameters 
% resultFile_wo_kanR_param_name_list = [resultFile_wo_kanR.estimated_param_summary.Name]; 
% resultFile_w_kanR_TX_param_name_list = [resultFile_w_kanR_TX.estimated_param_summary.Name];
% resultFile_w_kanR_TXTL_param_name_list = [resultFile_w_kanR_TXTL.estimated_param_summary.Name]; 
% 
% wo_kanR_param_summary = get_all_estimated_params(resultFile_wo_kanR.all_fitResults);
% w_kanR_TX_param_summary = get_all_estimated_params(resultFile_w_kanR_TX.all_fitResults);
% w_kanR_TXTL_param_summary = get_all_estimated_params(resultFile_w_kanR_TXTL.all_fitResults); 
% 
% all_param_bounds = [resultFile_wo_kanR.estimated_param_summary.Bounds];
% 
% figure; 
% for param_idx = 1:length(resultFile_wo_kanR_param_name_list)
% 
%     param_name = resultFile_wo_kanR_param_name_list{param_idx};
%     param_idx_kanR_TX = strcmp(resultFile_w_kanR_TX_param_name_list,param_name); 
%     param_idx_kanR_TXTL = strcmp(resultFile_w_kanR_TXTL_param_name_list,param_name); 
%     display_param_name = strrep(param_name,'_',' '); 
%     if length(display_param_name) > 15
%         display_param_name = sprintf('%s \n %s',display_param_name(1:15),display_param_name(16:end));
%     end
% 
%     % Plot param distribution for wo kanR
%     subplot(3,length(resultFile_wo_kanR_param_name_list),param_idx)
%     histogram(wo_kanR_param_summary(param_idx,:))
%     xlim(all_param_bounds(param_idx,:))
%     title(display_param_name)
%     if isequal(param_idx,1)
%         ylabel('w/o kanR')
%     end
% 
%     % Plot param distribution for w kanR TX 
%     subplot(3,length(resultFile_wo_kanR_param_name_list),length(resultFile_wo_kanR_param_name_list) * 1 + param_idx)
%     histogram(w_kanR_TX_param_summary(param_idx_kanR_TX,:))
%     xlim(all_param_bounds(param_idx,:))
%     if isequal(param_idx,1)
%         ylabel('w/ kanR TX')
%     end
% 
%     % Plot param distribution for w kanR TXTL 
%     subplot(3,length(resultFile_wo_kanR_param_name_list),length(resultFile_wo_kanR_param_name_list) * 2 + param_idx)
%     histogram(w_kanR_TXTL_param_summary(param_idx_kanR_TXTL,:))
%     xlim(all_param_bounds(param_idx,:))
%     if isequal(param_idx,1)
%         ylabel('w/ kanR TXTL')
%     end
% 
% end

% % Put all estimated parameters in the same place and plot SSE relationship
% % with individual parameter values in a scatter plot
% figure; 
% for param_idx = 1:length(resultFile_wo_kanR_param_name_list)
%     subplot(ceil(sqrt(length(resultFile_wo_kanR_param_name_list))),ceil(sqrt(length(resultFile_wo_kanR_param_name_list))),param_idx)
%     param_name = resultFile_wo_kanR_param_name_list{param_idx};
%     param_idx_kanR_TX = strcmp(resultFile_w_kanR_TX_param_name_list,param_name); 
%     param_idx_kanR_TXTL = strcmp(resultFile_w_kanR_TXTL_param_name_list,param_name); 
%     display_param_name = strrep(param_name,'_',' '); 
%     if length(display_param_name) > 15
%         display_param_name = sprintf('%s \n %s',display_param_name(1:15),display_param_name(16:end));
%     end
% 
%     scatter(wo_kanR_param_summary(param_idx,:),log10(resultFile_wo_kanR.metric_summary'),'black')
%     hold on
%     scatter(w_kanR_TX_param_summary(param_idx_kanR_TX,:),log10(resultFile_w_kanR_TX.metric_summary'),'blue')
%     scatter(w_kanR_TXTL_param_summary(param_idx_kanR_TXTL,:),log10(resultFile_w_kanR_TXTL.metric_summary'),'red')
%     title(display_param_name)
%     ylabel('log10(SSE)')
% 
% 
% end
% legend('without kanR','with kanR TX','with kanR TXTL')
%% Test posterior parameter sampling with PESTO toolbox 
%     % Use mRNA_prod_deg as test case
% load('param_est_run_save/20231019_mRNA_prod_deg_model_fitting_analysis.mat')
% 
%     % Construct parameter objects 
% param_name_list = {problemObject.Estimated.Name};
% all_param_min = nan(length(param_name_list),1); 
% all_param_max = nan(length(param_name_list),1); 
% all_param_guess = nan(length(param_name_list),1); 
% 
% [~,opt_idx] = min(metric_summary); 
% opt_fitResult = all_fitResults{opt_idx}; 
% for param_idx = 1:length(param_name_list)
%     param_obj = problemObject.Estimated(param_idx); 
%     all_param_min(param_idx) = param_obj.Bounds(1);
%     all_param_max(param_idx) = param_obj.Bounds(2);
%     all_param_guess(param_idx) = param_obj.InitialValue;
% end
% 
% parameters = struct('number',num2cell(length(param_name_list)),'min',all_param_min,'max',all_param_max,...
%     'name',{param_name_list'},'guess',all_param_guess);
% objective_Fun = @(paramsVec) temp_calculate_dev_from_data_mRNA_prod_deg(paramsVec); 
% 
%     % Construct objective function 
% [updated_parameters,fh] = getMultiStarts(parameters, objective_Fun); 



% [parameters,fh] = getParameterProfiles(parameters, objective_function, varargin); 
%% Screening for failure modes in large-scale fitting 
% data_table_file = load('param_est_run_save/20231005_large_scale_run_1002.mat'); 
% data_table = data_table_file.optimization_run_table; 
% 
    % Organize estimated parametres in one matrix and calculate parameter
    % correlation 
% all_estimated_params = nan(height(data_table),32); 
% for i = 1:height(data_table)
%     estimate_params = data_table{i,"ParameterEstimates"};
%     all_estimated_params(i,:) = estimate_params{1,1}; 
% end
%         % Load any problemObject to get the parameter names
% any_fitResult = load('param_est_run_save/20231102_param_est_run1534_92.mat');
% problemObject = any_fitResult.problemObject;
% param_names = {problemObject.Estimated.Name};
% param_names_for_lgd = strrep(param_names,'_',' '); 
% 
%         % Calculate correlation among estimated parameters 
% [corr_coef,corr_coef_p] = corrcoef(all_estimated_params); 
%             % Make nan diagonal values
% for param_idx = 1:size(corr_coef,1)
%     corr_coef(param_idx,param_idx) = nan;
% end
% 
%         % Use results from all fitting runs 
% figure; 
% h1 = heatmap(param_names_for_lgd,param_names_for_lgd,abs(corr_coef),'MissingDataColor','w');
% title('Parameter correlation in large-scale parameter fitting ')

        % Use filtered results 
% all_SSE = [data_table.SSE];
% median_SSE = median(all_SSE);
% keep_idx = find(all_SSE < median_SSE); 
% filtered_params = all_estimated_params(keep_idx,:);
% [filtered_corr_coef,filtered_corr_coef_p] = corrcoef(filtered_params); 
%             % Make nan diagonal values
% for param_idx = 1:size(filtered_corr_coef,1)
%     filtered_corr_coef(param_idx,param_idx) = nan;
% end
%         % Use results from above median runs only 
% figure; 
% h2 = heatmap(param_names_for_lgd,param_names_for_lgd,abs(filtered_corr_coef),'MissingDataColor','w');
% title('Parameter correlation in large-scale parameter fitting (with above median SSE)')

%     % Use top 10% SSE 
% sorted_SSE = sort(all_SSE,'ascend');
% SSE_threshold = sorted_SSE(480); 
% keep_idx = find(all_SSE < SSE_threshold); 
% filtered_params = all_estimated_params(keep_idx,:);
% [filtered_corr_coef,filtered_corr_coef_p] = corrcoef(filtered_params); 
%             % Make nan diagonal values
% for param_idx = 1:size(filtered_corr_coef,1)
%     filtered_corr_coef(param_idx,param_idx) = nan;
% end
%         % Use results from above median runs only 
% figure; 
% h3 = heatmap(param_names_for_lgd,param_names_for_lgd,abs(filtered_corr_coef),'MissingDataColor','w');
% title('Parameter correlation in large-scale parameter fitting (with top 10% SSE)')
% % 
% %     % Use cases that comply qualitatively 
% % filtered_params = all_estimated_params(exist_ideal_idx_list,:);
% % [filtered_corr_coef,filtered_corr_coef_p] = corrcoef(filtered_params); 
% %             % Make nan diagonal values
% % for param_idx = 1:size(filtered_corr_coef,1)
% %     filtered_corr_coef(param_idx,param_idx) = nan;
% % end
% %         % Use results from above median runs only 
% % figure; 
% % h3 = heatmap(param_names_for_lgd,param_names_for_lgd,abs(filtered_corr_coef),'MissingDataColor','w');
% % title('Parameter correlation in large-scale parameter fitting (qualitatively correct)')
% 
%     % Filter for cases where the positive crosstalk ratio for sigma70 weak
%     % is larger than the others 
% % Criteria: (1) The first ratio in sigma70 weak is larger than 1.1 (2) The
% % ratio is larger than the first ratio for other promotors 
% weak_prom_flag_array = nan(height(data_table),1); 
% for row_idx = 1:height(data_table)
%     crosstalk_ratios = data_table{row_idx,'CrosstalkRatios'}{1,1};
%     [T7_strong_crosstalk_ratio,T7_weak_crosstalk_ratio,sigma70_strong_crosstalk_ratio,sigma70_weak_crosstalk_ratio] = ...
%         unpack_crosstalk_ratios(crosstalk_ratios,promotor_name_list); 
%     weak_prom_flag_array(row_idx) = sigma70_weak_crosstalk_ratio(1) > 1.1 && ...
%  sigma70_weak_crosstalk_ratio(1) > T7_strong_crosstalk_ratio(1) + 0.1 && ...
%         sigma70_weak_crosstalk_ratio(1) > T7_weak_crosstalk_ratio(1) + 0.1 && ...
%                     sigma70_weak_crosstalk_ratio(1) > sigma70_strong_crosstalk_ratio(1) + 0.1;
% end

    % Correlation between SSE and sigma70_weak_avg+
% sigma70_weak_avg_pos = new_data_table.("sigma70_weak_+avg"); 
% filtered_sigma70_weak_avg_pos = sigma70_weak_avg_pos(~isnan(sigma70_weak_avg_pos)); 
% SSE = new_data_table.("SSE");
% filtered_SSE = SSE(~isnan(sigma70_weak_avg_pos)); 
% [rho,p_val] = corr(filtered_sigma70_weak_avg_pos,filtered_SSE); 

    % Is there any decrease in concentratio with really high plasmid
    % dosage? (This would involve running sims and would take long) 

% Load a problemObject
% any_fitResult = load('param_est_run_save/20231002_param_est_run1015_100.mat');
% problemObject = any_fitResult.problemObject;
% 
% % Create simFunction from the problemObject
% all_model_species = {problemObject.Model.Species.Name}; 
% reporter_species_name = 'protein sfGFP*'; 
% keep_idx = ~strcmp(all_model_species,reporter_species_name);
% additional_track_species = all_model_species(keep_idx); 
% data_table_species_name = [{reporter_species_name},additional_track_species];
% 
% simFunction = create_simFun_from_problemObject(problemObject,additional_track_species);
% dosing_information = create_dosing_info_from_problemObject(problemObject);
% 
% tStart = 0; 
% tEnd = 21600; 
% timepoint_oi = 10800; 

% % Run simFunction based on estimated parameters 
% flag_array = nan(height(data_table),4); 
% 
% for idx = 1:height(data_table)
%     estimated_param = data_table{idx,"ParameterEstimates"}{1,1};
%     label = data_table{idx,"Label"}{1,1}; 
%     [Time,Data] = simFunction(estimated_param',tEnd,dosing_information,tStart:tEnd);
%     % Criteria #1: Decrease in T7 strong, sigma70 strong, T7 weak baseline,
%     % no decrease in sigma70 weak baseline 
% 
%     % Get datapoints for the baselines at t=10800s. If there's a decrease
%     % with increasing plasmid concentrations, flag is true 
%     GFP_at_timepoint = nan(7,4); 
%     for prom_idx = 1:4
%         for data_idx = (prom_idx - 1) * 28 + 1:(prom_idx - 1) * 28 + 7 
%             selected_Time = Time{data_idx,1};
%             selected_Data = Data{data_idx,1}; 
%             GFP_val = get_GFP_at_timepoint(selected_Time,selected_Data,timepoint_oi);
%             GFP_at_timepoint(data_idx - 28 * (prom_idx - 1),prom_idx) = GFP_val; 
%         end
%     end
% 
%     GFP_at_timepoint_diff = diff(GFP_at_timepoint); 
%     flag_array(idx,:) = all(GFP_at_timepoint_diff > 0,1); 
% end

% load('param_est_run_save/20231011_large_scale_fitting_PEdecrease_check.mat'); 
% ideal_row = [0,0,0,1]; % Only sigma70_weak monotonically increases with plasmid dosage
% exist_ideal_loc = ones(size(flag_array,1),1); 
% for col_idx = 1:length(ideal_row)
%     ideal_num = ideal_row(col_idx); 
%     exist_ideal_loc = exist_ideal_loc & (flag_array(:,col_idx) == ideal_num); 
% end
% exist_ideal_idx_list = find(exist_ideal_loc); 
% row_oi = data_table(exist_ideal_idx_list,:); 
% 
% % Check if CUTP is the determining factor in those rows 
%     % Create simFunction
% all_species_name = {problemObject.Model.Species.Name}; 
% sfGFP_RNA_species_idx = contains(all_species_name,'RNA utrGFP--sfGFP') & ~contains(all_species_name,'RNase'); 
% sfGFP_RNA_species_track = all_species_name(sfGFP_RNA_species_idx); 
% resource_oi = [{'protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA'},sfGFP_RNA_species_track]; 
% resource_names = {'protein sfGFP','protein kanR','t7RNAP','RNAP','Ribo','RNase','AGTP','CUTP','AA','mRNA sfGFP'}; 
% reporter_species_name = 'protein sfGFP*'; 
% 
% simFunction = create_simFun_from_problemObject(problemObject,resource_oi);
% dosing_information = create_dosing_info_from_problemObject(problemObject);
% 
% tStart = 0; 
% tEnd = 21600; 
% timepoint_oi = 10800; 
% legend_names = {'0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM'};
% 
% for row_idx = 1:height(row_oi)
%     estimated_param = row_oi{row_idx,"ParameterEstimates"}{1,1}; 
%     [Time,Data] = simFunction(estimated_param',tEnd,dosing_information,tStart:tEnd);
% 
%     for prom_idx = 1:length(promotor_name_list)
%         promotor_name = promotor_name_list{prom_idx}; 
%         selected_Time = Time(28 * (prom_idx - 1) + 1:28 * (prom_idx - 1) + 7); 
%         selected_Data = Data(28 * (prom_idx - 1) + 1:28 * (prom_idx - 1) + 7); 
%         figure_title = sprintf('Large-scale fitting %d %s',row_idx,strrep(promotor_name,'_',' '));
%         save_path = sprintf('plots/20231016_large_scale_fitting_decreased_expression_resource_track_%s%s',row_oi{row_idx,"Label"},promotor_name);
%         plot_resource_time_course(selected_Time,selected_Data,resource_oi,resource_names,legend_names,figure_title,save_path)
%     end
% end

% %Looping through the data table to see whether any simulation captures
% %qualitatively accurate results 
%     % (1) At least one crosstalk ratio in the lowest 3 plasmid concentrations
%     % are larger than 1 
%     % (2) At least one crosstalk ratio in the highest 4 plasmid
%     % concentrations are smaller than 1
% flag_array_crosstalk_ratio = zeros(height(data_table),4); 
% for row_idx = 1:height(data_table)
%     crosstalk_ratios = data_table{row_idx,'CrosstalkRatios'}{1,1};
%     [T7_strong_crosstalk_ratio,T7_weak_crosstalk_ratio,sigma70_strong_crosstalk_ratio,sigma70_weak_crosstalk_ratio] = ...
%         unpack_crosstalk_ratios(crosstalk_ratios,promotor_name_list); 
%     for prom_idx = 1:length(promotor_name_list)
%         promotor_name = promotor_name_list{prom_idx};
%         eval(sprintf('crosstalk_ratio = %s_crosstalk_ratio;',promotor_name)); 
%         potential_pos_crosstalk_region = crosstalk_ratio(1:4,:); 
%         potential_neg_crosstalk_region = crosstalk_ratio(4:7,:); 
% 
%         if isequal(prom_idx,2) % In T7 weak there's no positive crosstalk observed in T7weak with T7 empty 
%             if any(any(potential_pos_crosstalk_region > 1.1,2)) && any(any(potential_neg_crosstalk_region < 0.9,2))
%                 flag_array_crosstalk_ratio(row_idx,prom_idx) = 1; 
%             end
%         else
%             if any(all(potential_pos_crosstalk_region > 1.1,2)) && any(all(potential_neg_crosstalk_region < 0.9,2))
%                 flag_array_crosstalk_ratio(row_idx,prom_idx) = 1; 
%             end
%         end
%     end
% end
% 
% ideal_row = [1,1,1,1]; % All promotors satisfy the condition 
% exist_ideal_loc = ones(size(flag_array_crosstalk_ratio,1),1); 
% for col_idx = 1:length(ideal_row)
%     ideal_num = ideal_row(col_idx); 
%     exist_ideal_loc = exist_ideal_loc & (flag_array_crosstalk_ratio(:,col_idx) == ideal_num); 
% end
% exist_ideal_idx_list = find(exist_ideal_loc); 
% row_oi = data_table(exist_ideal_idx_list,:); 

% For the crosstalk ratio I put, also plot out the time-course 
% load param_est_run_save/20231002_param_est_run1014_25.mat
% fitResult_oi = all_fitResults{26,1};
% simData = fitted(fitResult_oi);
% [sim_Time,sim_Data] = process_crosstalk_data_from_source(simData,'simulated_data'); 
% grouped_data = problemObject.Data;
% [exp_Time,exp_Data] = process_crosstalk_data_from_source(grouped_data,'grouped_data'); 
% plot_time_course_crosstalk(sim_Time,sim_Data,exp_Time,exp_Data)

% %Plot out crosstlak ratios for these 
% for oi_idx = 1:height(row_oi)
%     crosstalk_ratios = row_oi{oi_idx,"CrosstalkRatios"}{1,1}; 
%     label = row_oi{oi_idx,'Label'}; 
%     figure; 
%     for prom_idx = 1:length(promotor_name_list)
%         prom_crosstalk_ratio = crosstalk_ratios(:,prom_idx);
%         subplot(2,2,prom_idx)
%         scatter(conc_vec,prom_crosstalk_ratio{1},'MarkerFaceColor',[0.5,0.5,0.5])
%         hold on
%         scatter(conc_vec,prom_crosstalk_ratio{2},'MarkerFaceColor','r')
%         scatter(conc_vec,prom_crosstalk_ratio{3},'MarkerFaceColor','b')
% 
%         plot(conc_vec,prom_crosstalk_ratio{1},'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         plot(conc_vec,prom_crosstalk_ratio{2},'LineWidth',1.5,'Color','r')
%         plot(conc_vec,prom_crosstalk_ratio{3},'LineWidth',1.5,'Color','b')
%         title(strrep(promotor_name_list{prom_idx},'_',' '))
%     end
%     sgtitle(sprintf('Crosstalk ratio %s',strrep(label,'_',' ')))
%     legend('Empty','Empty T7','Empty sigma70')
%     saveas(gcf,sprintf('param_est_run_save/20231211_large_scale_crosstalk_ratio_%s',label),'png'); 
% end
% 
% % Check parameter distribution 
% all_params = nan(29,height(row_oi)); 
% for oi_idx = 1:height(row_oi)
%     estimated_params = row_oi{oi_idx,'ParameterEstimates'}{1,1}; 
%     all_params(:,oi_idx) = estimated_params; 
% end
% param_names = {problemObject.Estimated.Name};
% param_bounds = {problemObject.Estimated.Bounds}; 
% for param_idx = 1:length(param_names)
%     param_name = param_names{param_idx}; 
%     param_bound = param_bounds{param_idx}; 
%     subplot(5,6,param_idx)
%     histogram(all_params(param_idx,:))
%     hold on 
%     xline(param_bound(1),'LineWidth',1.5)
%     xline(param_bound(2),'LineWidth',1.5)
%     title(strrep(param_name,'_',' '))
% end
%% Quick characterization metrics for crosstalk ratios 
%     % Number of positive crosstalk ratios in T7strong, T7weak,
%     % sigma70strong, sigma70weak
%     % Average of positive crosstalk ratio/negative crosstalk ratio 
% 
% 
%     % Load the data table
%     data_table_file = load('param_est_run_save/20231005_large_scale_run_1002.mat'); 
%     data_table = data_table_file.optimization_run_table; 
% 
%     % Preassign a metric matrix 
%     metric_titles = {'SSE','LogLikelihood','total_+','total_-','total_+avg','total_-avg',...
%         'T7_strong_+','T7_weak_+','sigma70_strong_+','sigma70_weak_+',...
%         'T7_strong_-','T7_weak-','sigma70_strong_-','sigma70_weak_-',...
%         'T7_strong_+avg','T7_weak_+avg','sigma70_strong_+avg','sigma70_weak_+avg',...
%         'T7_strong_-avg','T7_weak_-avg','sigma70_strong_-avg','sigma70_weak_-avg'}; 
%     metric_mat = nan(height(data_table),length(metric_titles)); 
% 
%     % Calculate metrics based on crosstalk ratio cells 
%     for row_idx = 1:height(data_table)
%         crosstalk_ratio_cell = data_table.CrosstalkRatios(row_idx); 
%         crosstalk_ratio = crosstalk_ratio_cell{1,1}; 
%         SSE = data_table.SSE(row_idx); 
%         LogLikelihood = data_table.LogLikelihood(row_idx); 
%         % Fill the metric columns
%         metric_mat(row_idx,1) = SSE;
%         metric_mat(row_idx,2) = LogLikelihood; 
%         % Calculate crosstalk related metric 
%         prom_num_pos = zeros(1,4); 
%         prom_num_neg = zeros(1,4); 
%         prom_pos_crosstalk_avg = zeros(1,4);
%         prom_neg_crosstalk_avg = zeros(1,4); 
%         for prom_idx = 1:size(crosstalk_ratio,2)
%             crosstalk_ratio_prom = crosstalk_ratio(:,prom_idx); 
%             crosstalk_ratio_prom_mat = [crosstalk_ratio_prom{:}]; 
%             crosstalk_ratio_prom_mat_flattened = crosstalk_ratio_prom_mat(:); 
% 
%             pos_crosstalk_avg = mean(crosstalk_ratio_prom_mat_flattened(crosstalk_ratio_prom_mat_flattened > 1)); 
%             neg_crosstalk_avg = mean(crosstalk_ratio_prom_mat_flattened(crosstalk_ratio_prom_mat_flattened < 1)); 
% 
%             prom_num_pos(1,prom_idx) = sum(crosstalk_ratio_prom_mat_flattened > 1); 
%             prom_num_neg(1,prom_idx) = sum(crosstalk_ratio_prom_mat_flattened < 1); % ignoring when equals 1 - rarely happens due to numerical precision 
% 
%             prom_pos_crosstalk_avg(1,prom_idx) = pos_crosstalk_avg; 
%             prom_neg_crosstalk_avg(1,prom_idx) = neg_crosstalk_avg; 
%         end
%         total_num_pos = sum(prom_num_pos);
%         total_num_neg = sum(prom_num_neg); 
%         total_pos_avg = mean(prom_pos_crosstalk_avg,"all",'omitnan'); 
%         total_neg_avg = mean(prom_neg_crosstalk_avg,"all",'omitnan'); 
% 
%         % Assign to metric matrix 
%         metric_mat(row_idx,3:6) = [total_num_pos,total_num_neg,total_pos_avg,total_neg_avg];
%         metric_mat(row_idx,7:10) = prom_num_pos;
%         metric_mat(row_idx,11:14) = prom_num_neg;
%         metric_mat(row_idx,15:18) = prom_pos_crosstalk_avg;
%         metric_mat(row_idx,19:22) = prom_neg_crosstalk_avg;
%     end
% new_data_table = array2table(metric_mat,"VariableNames",metric_titles);
% new_data_table.Label = data_table.Label;
% 
%     % Apply criteria to 'filter' the new data table
% filtered_data_table_idx = new_data_table.('T7_strong_+') > 0 & new_data_table.('T7_weak_+') > 0 ...
%     & new_data_table.('sigma70_strong_+') > 0 & new_data_table.('sigma70_weak_+') > 0; 
% filtered_data_table = new_data_table(filtered_data_table_idx,:); 

%% Crosstalk ratio fitting for sigma70 weak 
% load param_est_run_save/20231004_param_est_run1246.mat
% resnorm_summary = nan(length(all_fitResults),1); 
% estimated_param_summary = nan(length(all_fitResults{1,1}.estimated_params),length(all_fitResults)); 
% for iter = 1:length(all_fitResults)
%     fitResult = all_fitResults{iter,1}; 
%     resnorm_summary(iter) = fitResult.resnorm; 
%     estimated_param_summary(:,iter) = fitResult.estimated_params; 
% end
% [min_resnorm,opt_idx] = min(resnorm_summary);
% opt_param = estimated_param_summary(:,opt_idx); 
% 
% % dosing info
% reporter_prom_name_list = {'DNA pJ23105--utrGFP--sfGFP'}; 
% conc_vec = [0.5,1,2.5,5,10,15,30]; 
% empty_conc = 10; 
% num_promotor = 1; 
% num_conc = 7; 
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
% % Play with parameters 
% param_names = [simFunction.Parameters.Name];
% tStart = 0;
% tEnd = 21600; 
% modify_param_name_list = {'TXTL_PJ23105_RNAPbound_R'}; 
% modify_param_val_list = [11.5078292516291]; 
% modified_params = modify_params(param_names,opt_params,modify_param_name_list,modify_param_val_list); 
% [Time,Data] = simFunction(modified_params',tEnd,dosing_information,tStart:tEnd);
% crosstalk_ratio = calculate_crosstalk_ratio_v2(Time,Data,num_conc,num_promotor,'PE');
% figure; 
% for conc_idx = 1:length(conc_vec)
%     subplot(3,3,conc_idx)
%     plot(Time{conc_idx,1},Data{conc_idx,1},'LineWidth',1.5,'Color','g')
%     hold on 
%     plot(Time{conc_idx + 1 * length(conc_vec),1},Data{conc_idx + 1 * length(conc_vec),1},'LineWidth',1.5,'Color',[0.5,0.5,0.5])
%     plot(Time{conc_idx + 2 * length(conc_vec),1},Data{conc_idx + 2 * length(conc_vec),1},'LineWidth',1.5,'Color','r')
%     plot(Time{conc_idx + 3 * length(conc_vec),1},Data{conc_idx + 3 * length(conc_vec),1},'LineWidth',1.5,'Color','b')
%     title(sprintf('%.1f nM',conc_vec(conc_idx)))
% end
% 
% 
% 
% % for iter = 1:48
% %     param = estimated_param_summary(:,iter); 
% %     [Time,Data] = simFunction(param',tEnd,dosing_information,tStart:tEnd);
% %     crosstalk_ratio = calculate_crosstalk_ratio_v2(Time,Data,num_conc,num_promotor,'PE');
% %     figure; 
% %     for conc_idx = 1:length(conc_vec)
% %         subplot(3,3,conc_idx)
% %         plot(Time{conc_idx,1},Data{conc_idx,1},'LineWidth',1.5,'Color','g')
% %         hold on 
% %         plot(Time{conc_idx + 1 * length(conc_vec),1},Data{conc_idx + 1 * length(conc_vec),1},'LineWidth',1.5,'Color',[0.5,0.5,0.5])
% %         plot(Time{conc_idx + 2 * length(conc_vec),1},Data{conc_idx + 2 * length(conc_vec),1},'LineWidth',1.5,'Color','r')
% %         plot(Time{conc_idx + 3 * length(conc_vec),1},Data{conc_idx + 3 * length(conc_vec),1},'LineWidth',1.5,'Color','b')
% %         title(sprintf('%.1f nM',conc_vec(conc_idx)))
% %     end
% %     sgtitle(sprintf('PE sigma70 weak crosstalk - ratio as obj fun (iter #%d)',iter))
% % end
%% Check rate mechanism fitting with tanh
% % load param_est_run_save/20231004_param_est_run1427_.mat
% % observable_name = 'protein sfGFP*'; 
% [~,opt_idx] = min(metric_summary);
% opt_fitResult = all_fitResults{opt_idx,1}; 
% eval(sprintf('opt_estimated_params = [estimated_param_summary.Estimate_%d];',opt_idx))
% estimated_param_names = [estimated_param_summary.Name]; 
% additional_track_species = {'AGTP','CUTP','AA','RNAP','t7RNAP','Ribo','toxin'}; 
% simFunction = create_simFun_from_problemObject(problemObject,additional_track_species);
% tStart = 0;
% tEnd = 21600; 
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% 
% % Modify parameters 
% modify_param_name_list = {'tx_capacity_param','AGTP_deg_time_0','toxin_threshold','k_toxin'};
% modify_param_value_list = [0.0001,4000,100,1000]; 
% modified_params = modify_params(estimated_param_names,opt_estimated_params,modify_param_name_list,modify_param_value_list); 
% 
% % Simulate 
% [Time,Data] = simFunction(modified_params',tEnd,dosing_information,tStart:tEnd);
% 
% % Plot out 
% figure; 
% conc_vec = [0.5,1,2.5,5,10,15,30]; 
% track_species_name = {'sfGFP','AGTP','CUTP','AA','RNAP','t7RNAP','Ribo','toxin'};  
% for conc_idx = 1:length(conc_vec)
%     Time_single = Time{conc_idx};
%     Data_single = Data{conc_idx}; 
%     for species_idx = 1:length(track_species_name)
%         species_name = track_species_name{species_idx}; 
%         subplot(3,3,species_idx)
%         plot(Time_single,Data_single(:,species_idx),'LineWidth',1.5)
%         hold on 
%         title(species_name)
% 
%     end
% end
% legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')

%% Play with tanh 
% figure; 
% x = -10:0.01:10; 
% y = 1 - tanh(10 * x); 
% plot(x,y)


% all_modes = {'TX','PE','combined'};
% for mode_idx = 1:length(all_modes)
%     mode = all_modes{mode_idx}; 
%     base_fitResult_file = load(sprintf('param_est_run_save/20230813_param_est_run1010_%s_fitResult1.mat',mode));
%     base_problemObject = base_fitResult_file.problemObject; 
%     for rep = 2:48
%         fitResult_file = load(sprintf('param_est_run_save/20230813_param_est_run1010_%s_fitResult%d.mat',mode,rep));
%         problemObject = fitResult_file.problemObject; 
%         % Check consistency for following: (1) length of data (2) error
%         % model (3) weights 
%         if ~strcmp(base_problemObject.ErrorModel,problemObject.ErrorModel)
%             fprintf('fit result not consistent in %s 1010',mode)
%             break
%         end
%         if ~all(isequal(problemObject.Weights,base_problemObject.Weights))
%             fprintf('fit result not consistent in %s 1010',mode)
%             break
%         end
%     end
% end