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
reporter_species_TX_name = 'RNA utrbroc--no_protein'; 

%% Alternative toxin mechanism result check 
% TX_gen_on_RNAP_resultFileName ='param_est_run_save/20231101_param_est_run0929_RNAP_toxin_mech.mat';
% TX_gen_on_Ribo_resultFileName = 'param_est_run_save/20230125_param_est_run1111_TXgen_Ribo.mat';
% TL_gen_on_RNAP_resultFileName = 'param_est_run_save/20230125_param_est_run2230_TLgen_RNAP.mat';
% TL_gen_on_Ribo_resultFileName = 'param_est_run_save/20230125_param_est_run1326_TLgen_Ribo.mat';
% 
% TX_gen_on_RNAP_resultFile = load(TX_gen_on_RNAP_resultFileName);
% TX_gen_on_Ribo_resultFile = load(TX_gen_on_Ribo_resultFileName); 
% TL_gen_on_RNAP_resultFile = load(TL_gen_on_RNAP_resultFileName);
% TL_gen_on_Ribo_resultFile = load(TL_gen_on_Ribo_resultFileName);
% 
% figure; 
% for iter = 1:length(TX_gen_on_RNAP_resultFile.all_fitResults)
%     TX_gen_on_RNAP_fitResult_oi = TX_gen_on_RNAP_resultFile.all_fitResults{iter};
%     fitted_data_TX_gen_on_RNAP = fitted(TX_gen_on_RNAP_fitResult_oi); 
%     [Time_TX_gen_on_RNAP,Data_TX_gen_on_RNAP] = process_crosstalk_data_from_source(fitted_data_TX_gen_on_RNAP,'simulated_data');
%     subplot(7,7,iter)
%     for conc_idx = 1:length(Data_TX_gen_on_RNAP)
%         plot(Time_TX_gen_on_RNAP{conc_idx}./3600,Data_TX_gen_on_RNAP{conc_idx},'LineWidth',1.5)
%         hold on 
%     end
% end
% sgtitle('TX on RNAP')
% 
% figure; 
% for iter = 1:length(TL_gen_on_Ribo_resultFile.all_fitResults)
%     TL_gen_on_Ribo_fitResult_oi = TL_gen_on_Ribo_resultFile.all_fitResults{iter};
%     fitted_data_TL_gen_on_Ribo = fitted(TL_gen_on_Ribo_fitResult_oi); 
%     [Time_TL_gen_on_Ribo,Data_TL_gen_on_Ribo] = process_crosstalk_data_from_source(fitted_data_TL_gen_on_Ribo,'simulated_data');
%     subplot(7,7,iter)
%     for conc_idx = 1:length(Data_TL_gen_on_Ribo)
%         plot(Time_TL_gen_on_Ribo{conc_idx}./3600,Data_TL_gen_on_Ribo{conc_idx},'LineWidth',1.5)
%         hold on 
%     end
% end
% sgtitle('TL on Ribo')
% 
% figure; 
% for iter = 1:length(TL_gen_on_RNAP_resultFile.all_fitResults)
%     TL_gen_on_RNAP_fitResult_oi = TL_gen_on_RNAP_resultFile.all_fitResults{iter};
%     fitted_data_TL_gen_on_RNAP = fitted(TL_gen_on_RNAP_fitResult_oi); 
%     [Time_TL_gen_on_RNAP,Data_TL_gen_on_RNAP] = process_crosstalk_data_from_source(fitted_data_TL_gen_on_RNAP,'simulated_data');
%     subplot(7,7,iter)
%     for conc_idx = 1:length(Data_TL_gen_on_RNAP)
%         plot(Time_TL_gen_on_RNAP{conc_idx}./3600,Data_TL_gen_on_RNAP{conc_idx},'LineWidth',1.5)
%         hold on 
%     end
% end
% sgtitle('TL on RNAP')
% 
% 
% figure; 
% for iter = 1:length(TX_gen_on_Ribo_resultFile.all_fitResults)
%     TX_gen_on_Ribo_fitResult_oi = TX_gen_on_Ribo_resultFile.all_fitResults{iter};
%     fitted_data_TX_gen_on_Ribo = fitted(TX_gen_on_Ribo_fitResult_oi); 
%     [Time_TX_gen_on_Ribo,Data_TX_gen_on_Ribo] = process_crosstalk_data_from_source(fitted_data_TX_gen_on_Ribo,'simulated_data');
%     subplot(7,7,iter)
%     for conc_idx = 1:length(Data_TX_gen_on_Ribo)
%         plot(Time_TX_gen_on_Ribo{conc_idx}./3600,Data_TX_gen_on_Ribo{conc_idx},'LineWidth',1.5)
%         hold on 
%     end
% end
% sgtitle('TX on Ribo')

%% Is the number of sample enough? 

%     % Load result file 
% result_file = load('test_save_files/202410_high_res_sampling_result_summary_updated.mat'); 
% all_penalty_terms = result_file.all_high_res_penalty_term; 
% valid_penalty_terms = all_penalty_terms(all(~isnan(all_penalty_terms),2),:); 
% init_param_file = load('test_save_files/20240819_lhs_sampled_params.mat'); 
% all_sampled_params = init_param_file.satisfied_sampled_params;
% valid_sampled_params = all_sampled_params(all(~isnan(all_penalty_terms),2),:);
% 
% N = 50000; 
% num_draft = 10; 
% 
% drafted_bs_sampled_params = cell(num_draft,1);
% drafted_bs_penalty_terms = cell(num_draft,1); 
% 
%     % Draft parameters and penalty terms 
% for draw_idx = 1:num_draft
%         % Draw N samples 
%     lottery = randperm(size(valid_penalty_terms,1)); 
%     select_idx_list = lottery(1:N); 
% 
%         % Get bootstrapped sampled params and penalty terms 
%     bs_sampled_params = valid_sampled_params(select_idx_list',:);
%     bs_penalty_terms = valid_penalty_terms(select_idx_list,:); 
% 
%         % Assign 
%     drafted_bs_sampled_params{draw_idx} = bs_sampled_params;
%     drafted_bs_penalty_terms{draw_idx} = bs_penalty_terms; 
% 
% end
% 
%     % Compare drafted penalty term pairwise 
%         % (No need to compare drafted parameters because we know they're
%         % from the same distribution) 
% num_penalty_term = size(bs_penalty_terms,2); 
% all_ks_test_h = cell(num_penalty_term,1);
% all_ks_test_p = cell(num_penalty_term,1); 
% 
% for penalty_idx = 1:num_penalty_term
% 
%     % initialize matrix to store num_draft * num_draft h's and p's 
%     penalty_ks_test_h = nan(num_draft,num_draft);
%     penalty_ks_test_p = nan(num_draft,num_draft); 
%     for draw_idx_1 = 1:num_draft
%         bs_penalty_terms_1 = drafted_bs_penalty_terms{draw_idx_1}; 
%         for draw_idx_2 = 1:num_draft
%             bs_penalty_terms_2 = drafted_bs_penalty_terms{draw_idx_2}; 
% 
%             bs_penalty_term_1 = bs_penalty_terms_1(:,penalty_idx);
%             bs_penalty_term_2 = bs_penalty_terms_2(:,penalty_idx); 
%             [h,p] = kstest2(bs_penalty_term_1,bs_penalty_term_2); 
%             penalty_ks_test_h(draw_idx_1,draw_idx_2) = h; 
%             penalty_ks_test_p(draw_idx_1,draw_idx_2) = p; 
%         end
%     end 
%     all_ks_test_h{penalty_idx} = penalty_ks_test_h; 
%     all_ks_test_p{penalty_idx} = penalty_ks_test_p; 
% 
%     % Plot out a heatmap for each 
%     figure;
%     heatmap(penalty_ks_test_h)
%     % sgtitle(strrep(result_file.penalty_labels{penalty_idx},'_',' '))
% 
% end

%% Figure S7
% result_file_name_wo_weights = 'param_est_run_save/20240829_param_est_run1207_wo_weight.mat';
% result_file = load(result_file_name_wo_weights); 
% dosing_information = create_dosing_info_from_problemObject(result_file.problemObject);
% % 
% % [~,opt_idx] = max(result_file.all_SSE); 
% % opt_fitResult = result_file.all_fitResults{opt_idx}; 
% % fitted_data = fitted(opt_fitResult);
% % [simulated_time,simulated_data] = process_crosstalk_data_from_source(fitted_data,'simulated_data');
% % plot_simulated(simulated_time,simulated_data,'crosstalk_ratio'); 
% % plot_simulated(simulated_time,simulated_data,'baseline'); 
% 
% 
% 
% 
% all_fitResults = result_file.all_fitResults; 
% 
% all_high_res_penalty_terms = nan(length(all_fitResults),32);
% all_SSE = nan(length(all_fitResults),1);
% for iter = 1:length(all_fitResults)
% 
%     fitResult = all_fitResults{iter}; 
%     [scaled_high_res_penalty_terms,~,~,~] = wrapper_calculate_obj_higher_resolution([fitResult.ParameterEstimates.Estimate]',...
%         result_file.problemObject,dosing_information,[],[],mode,num_conc,num_promotor); 
%     all_high_res_penalty_terms(iter,:) = scaled_high_res_penalty_terms; 
%     all_SSE(iter) = fitResult.SSE; 
% 
% 
% end
% 
% save(result_file_name_wo_weights,'all_high_res_penalty_terms','all_SSE','-append')

result_file_name_w_weights = 'param_est_run_save/20240829_param_est_run1207_w_weight.mat';
result_file = load(result_file_name_w_weights); 
dosing_information = create_dosing_info_from_problemObject(result_file.problemObject);


all_fitResults = result_file.all_fitResults; 

all_high_res_penalty_terms = nan(length(all_fitResults),32);
all_SSE = nan(length(all_fitResults),1);
for iter = 1:length(all_fitResults)

    fitResult = all_fitResults{iter}; 
    [scaled_high_res_penalty_terms,~,~,~] = wrapper_calculate_obj_higher_resolution([fitResult.ParameterEstimates.Estimate]',...
        result_file.problemObject,dosing_information,[],[],mode,num_conc,num_promotor); 
    all_high_res_penalty_terms(iter,:) = scaled_high_res_penalty_terms; 
    all_SSE(iter) = fitResult.SSE; 


end

save(result_file_name_w_weights,'all_high_res_penalty_terms','all_SSE','-append')
% 
% all_SSE = result_file.all_SSE;
% all_high_res_penalty_terms = result_file.all_high_res_penalty_terms; 
% 
% high_res_penalty_goals = load('test_save_files/high_resolution_penalty_goals.mat'); 
% penalty_cutoff = high_res_penalty_goals.penalty_goals; 


% % Visualize the SSE vs. crosstalk ratio penalties 
% % [sorted_SSE,sort_idx] = sort(all_SSE,'ascend');
% % sorted_all_high_res_penalty_terms = all_high_res_penalty_terms(sort_idx,:); 
% % positive_crosstalk_penalty = sorted_all_high_res_penalty_terms(:,7:18);
% % negative_crosstalk_penalty = sorted_all_high_res_penalty_terms(:,19:30); 
% % figure; 
% % yline(-0.05,'LineWidth',1.5)
% % hold on 
% all_num_satisfied_penalty = nan(length(all_SSE),1); 
% for iter = 1:length(all_SSE)
%     % x_SSE = log10(sorted_SSE(iter)); 
%     % selected_positive_crosstalk_penalty = positive_crosstalk_penalty(iter,:);
%     % selected_negative_crosstalk_penalty = negative_crosstalk_penalty(iter,:);
%     % positive_satisfied = selected_positive_crosstalk_penalty(selected_positive_crosstalk_penalty <= -0.05); 
%     % positive_not_satisfied = selected_positive_crosstalk_penalty(selected_positive_crosstalk_penalty > -0.05); 
%     % negative_satisfied = selected_negative_crosstalk_penalty(selected_negative_crosstalk_penalty <= -0.05); 
%     % negative_not_satisfied = selected_negative_crosstalk_penalty(selected_negative_crosstalk_penalty > -0.05); 
% 
%     num_satisfied_penalty = sum(sorted_all_high_res_penalty_terms(iter,:) < penalty_cutoff); 
%     all_num_satisfied_penalty(iter) = num_satisfied_penalty; 
% 
%     % scatter(x_SSE,length(positive_satisfied) + length(negative_satisfied),'o','filled')
%     % hold on 
%     % scatter(x_SSE,length(positive_not_satisfied),'o')
%     % scatter(x_SSE,length(negative_satisfied),"^",'filled')
%     % scatter(x_SSE,length(negative_not_satisfied),"^")
% 
% end
% scatter(log10(sorted_SSE),all_num_satisfied_penalty,'filled'); 
% hold on 
% all_SSE = nan(96,1);
% % Check initial points 
% problemObject = result_file.problemObject; 
% all_initial_estimated_params = nan(length(all_SSE),31); 
% all_init_penalty_term = nan(length(all_SSE),32); 
% all_init_SSE = nan(length(all_SSE),1); 
% 
% for iter = 1:length(all_SSE)
% 
%     rng(iter)
%     % Use initial parameter values for first run 
%     if ~isequal(iter,1)
%         % Then sample initial parameter values within parameter bounds
%         % in the log space 
%         InitialValue = nan(1,31); 
%         for est_param_idx = 1:length(problemObject.Estimated)
%             InitialTransformedValue =  problemObject.Estimated(est_param_idx).TransformedBounds(1) + ...
%                 (problemObject.Estimated(est_param_idx).TransformedBounds(2) -  problemObject.Estimated(est_param_idx).TransformedBounds(1)) * rand;
%             InitialValue(est_param_idx) =  exp(InitialTransformedValue);
%         end
%     else
%         InitialValue = [problemObject.Estimated.InitialValue]; 
%     end
%     all_initial_estimated_params(iter,:) = InitialValue; 
% 
%     [scaled_high_res_penalty_terms,~,~,~,SSE] = wrapper_calculate_obj_higher_resolution(InitialValue,...
%     problemObject,dosing_information,[],[],mode,num_conc,num_promotor); 
%     all_init_penalty_term(iter,:) = scaled_high_res_penalty_terms; 
%     all_init_SSE(iter) = SSE; 
% 
% end
% save(result_file_name_wo_weights,'all_init_penalty_term','all_init_SSE','-append')

% scatter(log10(all_init_SSE),sum(all_init_penalty_term < penalty_cutoff,2))

%% Update the problemObject in previous init param file 
% init_param_file = load('test_save_files/20240819_lhs_sampled_params.mat'); 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2_expanded.xlsx';
% problemObject_updated = setProblemObject_v3_OctUpdate(group_description,[],param_info_path); 
% save('test_save_files/20240819_lhs_sampled_params.mat','problemObject_updated','-append')

%% Test updated mechanism with macromolecule toxicity 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2_expanded.xlsx';
% problemObject = setProblemObject_v3_OctUpdate(group_description,[],param_info_path); 
% save_path_id = 'macromol_toxin';
% wrapper_fit_to_data_w_probObject_v2(problemObject,save_path_id)

% result_file_macromol_toxin = load('param_est_run_save/20241028_param_est_run1346_macromol_toxin.mat'); 
% all_fitResult = result_file_macromol_toxin.all_fitResults; 
% 
% all_SSE = nan(length(all_fitResult),1);
% all_estimated_params = nan(length(all_fitResult),height(all_fitResult{1}.ParameterEstimates)); 
% for iter = 1:length(all_fitResult)
%     fitResult = all_fitResult{iter};
%     if ~isempty(fitResult)
%         all_SSE(iter) = fitResult.SSE; 
%         all_estimated_params(iter,:) = fitResult.ParameterEstimates.Estimate'; 
%     end
% end

    % [~,min_idx] = min(all_SSE); 
    % fitResult_oi = all_fitResult{min_idx}; 

% figure;
% for iter = 1:100
%     subplot(10,10,iter)
%     fitResult_oi = all_fitResult{iter}; 
%     if ~isempty(fitResult_oi)
%         simData_oi = fitted(fitResult_oi); 
%         for data_idx = 1:length(simData_oi)
%             simData_oi_single = simData_oi(data_idx); 
%             sfGFP_name_idx = strcmp(simData_oi_single.DataNames,reporter_species_name);
%             sim_time = simData_oi_single.Time; 
%             sim_data = simData_oi_single.Data(:,sfGFP_name_idx); 
%             plot(sim_time,sim_data,'LineWidth',1.5)
%             hold on 
% 
%         end
%     end
% end

% % Close study on iterations where the decrease is captured 
% idx_oi_list = [8,22,32,35,44]; 
% for idx_oi = idx_oi_list
%     fitResult_oi = all_fitResult{idx_oi}; 
%     fitted_data = fitted(fitResult_oi); 
% 
%     % Plot out all species for 15 nM & 30 nM 
%     simData_oi_15nM = fitted_data(6);
%     simData_oi_30nM = fitted_data(7); 
%     figure; 
%     num_species = length(simData_oi_30nM.DataNames);
%     for species_idx = 1:num_species
%         subplot(ceil(sqrt(num_species)),ceil(sqrt(num_species)),species_idx)
%         plot(simData_oi_15nM.Time,simData_oi_15nM.Data(:,species_idx),'LineWidth',1.5)
%         hold on 
%         plot(simData_oi_30nM.Time,simData_oi_30nM.Data(:,species_idx),'LineWidth',1.5)
%         title(strrep(simData_oi_15nM.DataNames{species_idx},'_',' '))
%     end
% 
% 
% end

%% Flip positive & negative crosstalk ratio 
% load('test_save_files/20241023_model_oi_info.mat');
% % param_name_list = {problemObject.Estimated.Name}'; 
% [simulated_time,simulated_data] = simFunction(sampled_params_oi,tEnd,dosing_information,tStart:tEnd);
% species_name_list = simFunction.Observables.Name;
% species_option_struct = struct('promotor_oi','T7_strong','conc_idx',4);
% species_option_struct.species_name_list = species_name_list; 
% 
% plot_simulated(simulated_time,simulated_data,'BundledSpecies',[],[],species_option_struct); 
% plot_simulated(simulated_time,simulated_data,'crosstalk_ratio'); 
% plot_simulated(simulated_time,simulated_data,'baseline'); 

% TX_group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% param_info_path = 'parameters_v2_expanded.xlsx';
% problemObject_TX = setProblemObject_v3(TX_group_description,[],param_info_path);
% species_name_list = {problemObject_TX.Model.Species.Name};
% reporter_species_idx = find(strcmp(species_name_list,reporter_species_TX_name));
% additional_track_species_list = [species_name_list(1:reporter_species_idx - 1),species_name_list(reporter_species_idx + 1:end)];
% simFunction_TX = create_simFun_from_problemObject(problemObject_TX,additional_track_species_list);
% dosing_information_TX = create_dosing_info_from_problemObject(problemObject_TX);
% [simulated_time_TX,simulated_data_TX] = simFunction_TX(sampled_params_oi,tEnd,dosing_information_TX,tStart:tEnd);
% plot_simulated(simulated_time_TX,simulated_data_TX,'exhaustive'); 
% plot_time_course_crosstalk(simulated_time_TX,simulated_data_TX)
%% OAT sensitivity analysis on parameter set of interest 
% load('test_save_files/20241023_model_oi_info.mat'); 
% param_name_list = {problemObject.Estimated.Name}'; 
% param_bounds = [problemObject.Estimated.Bounds]; 
% param_bounds = reshape(param_bounds,[2,31]);
% param_bounds = param_bounds';
% log_param_bounds = log10(param_bounds); 
% log_sampled_param_oi = log10(sampled_params_oi); 
% 
% % lhs sampling
% num_sample = 100; 
% lhsSamples = lhsdesign(num_sample,1); 
% 
% all_high_res_penalty_terms_OATSA = nan(length(param_name_list) * num_sample,32); 
% for param_idx = 1:length(param_name_list)
%     log_param_bounds_single = log_param_bounds(param_idx,:); 
%     sampled_param_single = log_param_bounds_single(1) + (log_param_bounds_single(2) - log_param_bounds_single(1)) * lhsSamples;
%     nominal_sampled_param = 10 .^ sampled_param_single;
%     all_high_res_penalty_terms_OATSA_param = nan(num_sample,32); 
%     parfor sample_idx = 1:length(nominal_sampled_param)
%         modified_params = sampled_params_oi; 
%         modified_params(param_idx) = nominal_sampled_param(sample_idx);
%         [simulated_time,simulated_data] = simFunction(modified_params,tEnd,dosing_information,tStart:tEnd);
%         [scaled_penalty_vec,~,~,~] = wrapper_calculate_obj_higher_resolution_from_data(simulated_time,simulated_data,...
%                                                 problemObject,mode,num_conc,num_promotor);
%         all_high_res_penalty_terms_OATSA_param(sample_idx,:) = scaled_penalty_vec; 
% 
%     end
%     all_high_res_penalty_terms_OATSA((param_idx - 1) * num_sample + 1:param_idx * num_sample,:) = all_high_res_penalty_terms_OATSA_param; 
% end
% save('test_save_files/20241023_OAT_sensitivity_analysis_onParamSetoi.mat'); 
   
% load('test_save_files/20241023_OAT_sensitivity_analysis_onParamSetoi.mat'); 
% [baseline_high_res_penalty,~,~,~] = wrapper_calculate_obj_higher_resolution(sampled_params_oi,problemObject,dosing_information,[],[],mode,num_conc,num_promotor); 
% penalty_idx_list = [1:6,31:32];
% for param_idx = 1:length(param_name_list)
%     subplot(6,6,param_idx)
%     log_param_bounds_single = log_param_bounds(param_idx,:); 
%     sampled_param_single = log_param_bounds_single(1) + (log_param_bounds_single(2) - log_param_bounds_single(1)) * lhsSamples;
% 
%     for penalty_idx = penalty_idx_list%size(all_high_res_penalty_terms_OATSA,2)
%         baseline_penalty = baseline_high_res_penalty(penalty_idx);
%         vec_penalty_oi = all_high_res_penalty_terms_OATSA((param_idx - 1) * 100 + 1:param_idx * 100,penalty_idx); 
% 
%         [sorted_sampled_params,sort_idx] = sort(sampled_param_single,'ascend');
%         sorted_vec_penalty_oi = vec_penalty_oi(sort_idx); 
% 
%         plot(sorted_sampled_params,sorted_vec_penalty_oi)
%         hold on 
%         scatter(log10(sampled_params_oi(param_idx)),baseline_high_res_penalty(penalty_idx));
% 
% 
%     end
%     title(strrep(param_name_list{param_idx},'_',' '))
% end
% legend('baseline dev','crosstalk ratio dev','T7 strong baseline penal','T7 weak baseline penal','sigma70 strong baseline penal',...
%     'sigma70 weak baseline penal','Large crosstalk ratio','Residual mRNA')
%% Run missing parameter sampling results 

% parfor iter = 1:1591
%     result_file_name = sprintf('param_est_run_save/20241020_param_sampling_constrained_run%d.mat',run_idx); 
%     if ~exist(result_file_name,'file')
%         run_param_and_init_cond_sampling_v2(iter); 
%     end
% end

%% Fit TX-level crosstalk after updating empty degradation parameter

% TX_group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% param_info_path = 'parameters_v2_expanded.xlsx';
% problemObject_TX = setProblemObject_v3(TX_group_description,[],param_info_path);
% 
%     % Modify parameter bounds for degradation-related parameters? 
% RNAdeg_fitting_result = load('test_save_files/20240916_alternative_kinetics_fitting_original_data_truobleshooted_mechanism_part2.mat'); 
% 
% all_fval = RNAdeg_fitting_result.all_fval_MassAction; 
% filter_idx_list = find(all_fval < 1e+08); 
% all_estimated_params = RNAdeg_fitting_result.all_estimated_params_MassAction; 
% filtered_estimated_params = all_estimated_params(filter_idx_list,:); 

% 
% figure; 
% param_bounds_CI = nan(size(filtered_estimated_params,2),2); 
% for param_idx = 1:size(filtered_estimated_params,2)
%     test_log_params = log10(filtered_estimated_params(:,param_idx)); 
%     pd = fitdist(test_log_params,'Kernel'); 
%     pd_normal = fitdist(test_log_params,'normal'); 
% 
%     test_x = min(test_log_params) -1 :0.01:max(test_log_params) + 1;
%     pdf_normal = pdf(pd_normal,test_x);
%     pdf_kernel = pdf(pd,test_x);
%     subplot(2,2,param_idx)
%     plot(test_x,pdf_normal)
%     hold on 
%     plot(test_x,pdf_kernel)
% 
%     ci = paramci(pdf_normal);
%     param_bounds_CI(param_idx,:) = ci(:,1)'; 
% 
%     xline(ci(1,1))
%     xline(ci(2,1))
% end
% 
% 
% dosing_information_TX = create_dosing_info_from_problemObject(problemObject_TX); 

%% Organize & Analyze results 
%     % Organize result files into summary 
% result_fileName_example = 'param_est_run_save/20241029_param_sampling_constrained_run856.mat';
% result_file_example = load(result_fileName_example);
% 
% init_param_file = load('test_save_files/20240819_lhs_sampled_params.mat'); 
% num_result_file = ceil(size(init_param_file.satisfied_sampled_params,1) / 100); 
% all_sampled_params = init_param_file.satisfied_sampled_params; 
% 
% % Preassign space for penalty terms 
% all_high_res_penalty_term = nan(size(init_param_file.satisfied_sampled_params,1),size(result_file_example.high_res_all_penalty_terms,2));
% for iter = 1:num_result_file
% 
%     % Load result file
%     result_file_name = sprintf('param_est_run_save/20241029_param_sampling_constrained_run%d.mat',iter); 
%     if exist(result_file_name,'file')
%         result_file = load(result_file_name); 
% 
%         % Extract information 
%         sampled_params_selected = result_file.sampled_params_selected;
%         sampled_penalty = result_file.high_res_all_penalty_terms; 
%         % if ~exist('problemObject','var')
%         %     % problemObject = result_file.problemObject;
%         %     dosing_information = result_file.dosing_information;
%         %     penalty_term_labels = result_file.penalty_term_labels;
%         %     penalty_term_length = result_file.penalty_term_length; 
%         % end
% 
%         % Sanity check for sampled parameters alignment 
%         init_param_file_sampled_params = init_param_file.satisfied_sampled_params((iter - 1) * 100 + 1:iter * 100,:); 
%         if ~all(init_param_file_sampled_params == sampled_params_selected)
%             error('Mismatched sampled parameters')
%         end
% 
%         % Compile 
%         all_high_res_penalty_term((iter - 1) * 100 + 1:iter * 100,:) = sampled_penalty; 
%     end
% 
% end
% 
% save('test_save_files/202411_high_res_sampling_result_summary_updated.mat','all_high_res_penalty_term');

% result_summary_file = load('test_save_files/202410_high_res_sampling_result_summary_updated.mat'); 
% init_param_file = load('test_save_files/20240819_lhs_sampled_params.mat'); 
% sampled_result_file = load('param_est_run_save/20241020_param_sampling_constrained_run856.mat');
% 
%     % Load penalty cutoff 
% high_res_penalty_goals = load('test_save_files/high_resolution_penalty_goals.mat'); 
% penalty_cutoff = high_res_penalty_goals.penalty_goals; 
% 
%     % Create simFunction
% all_model_species = {init_param_file.problemObject.Model.Species.Name}; 
% reporter_species_idx = find(strcmp(all_model_species,reporter_species_name)); 
% additional_track_species = [all_model_species(1:reporter_species_idx - 1),all_model_species(reporter_species_idx + 1:end)];
% simFunction = create_simFun_from_problemObject(init_param_file.problemObject,additional_track_species); 
% 
% 
% satisfy_flag_matrix = streamline_penalty_criteria(result_summary_file.all_high_res_penalty_term,penalty_cutoff);
% 
% num_satisfied_criteria_per_sample = sum(satisfy_flag_matrix,2); 
% num_satisfied_sample_all = nan(size(satisfy_flag_matrix,2),1);
% for num_satisfied = 1:size(satisfy_flag_matrix,2)
%     num_satisfied_sample = sum(num_satisfied_criteria_per_sample >= num_satisfied); 
%     num_satisfied_sample_all(num_satisfied) = num_satisfied_sample; 
% end
% 
% funnelchart(num_satisfied_sample_all)

% num_satisfied_8_sample_idx_list = find(num_satisfied_criteria_per_sample >= 8); 
% high_res_penalty_terms_oi = result_summary_file.all_high_res_penalty_term(num_satisfied_8_sample_idx_list,:); 

% Find sampled parameter sets that satisfy key criteria
% sample_idx_list_oi = find(result_summary_file.all_high_res_penalty_term(:,1) < 11.4 & all(result_summary_file.all_high_res_penalty_term(:,3:30) < penalty_cutoff(3:30),2));
% 
% species_name_list = simFunction.Observables.Name;
% species_option_struct = struct('promotor_oi','sigma70_strong','conc_idx',7);
% species_option_struct.species_name_list = species_name_list; 
% 
% for k = 1:length(sample_idx_list_oi)
%     idx_oi = sample_idx_list_oi(k);
%     sampled_params_oi = init_param_file.satisfied_sampled_params(idx_oi,:); 
    % [simulated_time,simulated_data] = simFunction(sampled_params_oi,tEnd,init_param_file.dosing_information,tStart:tEnd);
% 
%     % plot_simulated(simulated_time,simulated_data,'baseline')
%     % plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
    % plot_simulated(simulated_time,simulated_data,'BundledSpecies',[],[],species_option_struct)
%     plot_simulated(simulated_time,simulated_data,'exhaustive')
% end
% 
% plot_qual_obj_dist(high_res_penalty_terms_oi,penalty_cutoff)

% TX_fitted_result_file_updated = load('param_est_run_save/20241008_param_est_run1356_TX.mat');
% 
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
%  % Plot out predicted time-course for every iteration
%     if iter <= 10
%         simData_temp = fitted(fitResult); 
%         [simulated_time_updated_temp,simulated_data_updated_temp] = ...
%             process_crosstalk_data_from_source(simData_temp,'simulated_data');
%         plot_simulated(simulated_time_updated_temp,simulated_data_updated_temp,'exhaustive')
%         % plot_time_course_crosstalk(simulated_time_updated_temp,simulated_data_updated_temp)   
%     end
% 
% end
% 
% 
% % TX-level integration 
% TX_group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% param_info_path = 'parameters_v2_expanded.xlsx';
% problemObject_TX = setProblemObject_v3(TX_group_description,[],param_info_path);
% dosing_information_TX = create_dosing_info_from_problemObject(problemObject_TX); 

%     % Create simFunction for TX
% all_model_species_TX = {problemObject_TX.Model.Species.Name}; 
% reporter_species_idx_TX = find(strcmp(all_model_species_TX,'RNA utrbroc--no_protein')); 
% additional_track_species_TX = [all_model_species_TX(1:reporter_species_idx_TX - 1),all_model_species_TX(reporter_species_idx_TX + 1:end)];
% simFunction_TX = create_simFun_from_problemObject(problemObject_TX,additional_track_species_TX); 
%     % Organize sampled parameters in protein-level crosstalk
% protein_crosstalk_param_names = simFunction.Parameters.Name; 
% TX_crosstalk_param_names = simFunction_TX.Parameters.Name;
% sampled_params_oi_TX = nan(length(TX_crosstalk_param_names),1); 
% for TX_param_idx = 1:length(sampled_params_oi_TX)
%     TX_param_name = TX_crosstalk_param_names{TX_param_idx};
%     param_idx = find(strcmp(protein_crosstalk_param_names,TX_param_name));
%     if ~isempty(param_idx)
%         sampled_params_oi_TX(TX_param_idx) = sampled_params_oi(param_idx); 
%     end
% 
% end
% 
% sampled_params_oi_TX(12) = 0.02; 
% [simulated_time,simulated_data] = simFunction_TX(sampled_params_oi_TX',tEnd,dosing_information_TX,tStart:tEnd);
% plot_simulated(simulated_time,simulated_data,'exhaustive')

%% Update TXTL_RNA_deg_R for empty vectors (should be same as sfGFP) 

% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2_expanded.xlsx';
% problemObject = setProblemObject_v3_OctUpdate(group_description,[],param_info_path); 

%% Analyze results from fitting SSE and qual obj 
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
% %     % Compile results for results using qualitative features as objective function
% %         % Record: init & final params + penalty, fval, 
% % num_params = length(fitting_result_file_SSE_sample.all_fitResults{1,1}.estimated_params); 
% % all_init_params = nan(16 * 100,num_params);
% % all_estimated_params = nan(16 * 100,num_params);
% % all_optimal_obj = nan(16 * 100,1); 
% % all_init_penalty = nan(16 * 100,7);
% % all_final_penalty = nan(16 * 100,7); 
% % all_exitFlag = nan(16 * 100,1); 
% % 
% % for idx = 1:16
% %         % Saving names were mixed up - this is counterintuitive but correct
% %     % fitting_result_file_name_qual_obj = sprintf('param_est_run_save/20240822_lhs_sampled_params_opt_run%d.mat',idx); 
% %     fitting_result_file_name_qual_obj = sprintf('param_est_run_save/20240822_lhs_sampled_params_opt_qual_obj_run%d.mat',idx);
% %     fitting_result_file_qual_obj= load(fitting_result_file_name_qual_obj);
% %     % fitting_result_file_qual_obj = load(fitting_result_file_name_qual_obj); 
% % 
% %     % % Sanity check with initial parameter file 
% %     % sampled_params_from_init_param_file = init_params_file_all_sampled_params((idx - 1) * 100 + 1:idx * 100,:);
% %     % sampled_params_from_result_file = fitting_result_file_qual_obj.; 
% %     % if ~all(all(sampled_params_from_init_param_file == sampled_params_from_result_file))
% %     %     error('Sanity check on initial parameter failed')
% %     % end
% % 
% %     % Access fields in result file 
% %     all_fitResults = fitting_result_file_qual_obj.all_fitResults;
% %     estimated_params_org = nan(length(all_fitResults),num_params);
% %     fval_org = nan(length(all_fitResults),1); 
% %     exitflag_org = nan(length(all_fitResults),1); 
% % 
% %     for iter = 1:length(all_fitResults)
% %         if ~isempty(all_fitResults{iter,1})
% %             estimated_params_org(iter,:) = all_fitResults{iter,1}.ParameterEstimates.Estimate'; 
% %             fval_org(iter) = all_fitResults{iter,1}.SSE;
% %             exitflag_org(iter) = all_fitResults{iter,1}.ExitFlag; 
% %         end
% %     end
% % 
% %     % Assign values 
% %     all_init_params((idx - 1) * 100 + 1:idx * 100,:) = init_params_file_all_sampled_params((idx - 1) * 100 + 1:idx * 100,:); 
% %     all_estimated_params((idx - 1) * 100 + 1:idx * 100,:) = estimated_params_org;
% %     all_optimal_obj((idx - 1) * 100 + 1:idx * 100,:) = fval_org; 
% %     % all_init_penalty((idx - 1) * 100 + 1:idx * 100,:) = fitting_result_file_qual_obj.all_penalty_term_init; 
% %     % all_final_penalty((idx - 1) * 100 + 1:idx * 100,:) = fitting_result_file_qual_obj.all_penalty_term_final; 
% %     all_exitFlag((idx - 1) * 100 + 1:idx * 100,:) = exitflag_org;
% % 
% % end
% % save('test_save_files/20240919_lhd_sampled_params_optimizatio_SSE_result_summary.mat','all_init_params','all_estimated_params','all_init_penalty',...
% %     'all_final_penalty','all_exitFlag');
% 
% result_file_SSE = load('test_save_files/20240919_lhd_sampled_params_optimizatio_SSE_result_summary.mat'); 
% result_file_qual_obj = load('test_save_files/20240919_lhd_sampled_params_optimizatio_qual_obj_result_summary.mat'); 
% 
% % % Calculate high resolution penalty terms
% 
% %     % Create simFunction and dosing information 
% simFunction = fitting_result_file_SSE_sample.simFunction; 
% problemObject = fitting_result_file_SSE_sample.problemObject; 
% dosing_information = create_dosing_info_from_problemObject(problemObject);
% % 
% % all_high_res_penalty_SSE = nan(size(result_file_SSE.all_estimated_params,1),32); 
% % all_high_res_penalty_qual_obj = nan(size(result_file_qual_obj.all_estimated_params,1),32); 
% % parfor iter = 1:size(result_file_SSE.all_estimated_params,1)
% %     estimated_params_oi_SSE = result_file_SSE.all_estimated_params(iter,:); 
% %     if ~all(isnan(estimated_params_oi_SSE))
% %         [scaled_penalty_vec_SSE,penalty_term_labels,penalty_term_length,~] = ...
% %             wrapper_calculate_obj_higher_resolution(estimated_params_oi_SSE,problemObject,dosing_information,[],[],mode,num_conc,num_promotor); 
% %         all_high_res_penalty_SSE(iter,:) = scaled_penalty_vec_SSE; 
% %     end
% % 
% %     estimated_params_oi_qual_obj = result_file_qual_obj.all_estimated_params(iter,:); 
% %     if ~all(isnan(estimated_params_oi_qual_obj))
% %         [scaled_penalty_vec_qual_obj,penalty_term_labels,penalty_term_length,~] = ...
% %             wrapper_calculate_obj_higher_resolution(estimated_params_oi_qual_obj,problemObject,dosing_information,[],[],mode,num_conc,num_promotor); 
% %         all_high_res_penalty_qual_obj(iter,:) = scaled_penalty_vec_qual_obj; 
% %     end
% % 
% % end
% % 
% % save('test_save_files/20240919_lhd_sampled_params_optimizatio_SSE_result_summary.mat','all_high_res_penalty_SSE','-append');
% % save('test_save_files/20240919_lhd_sampled_params_optimizatio_qual_obj_result_summary.mat','all_high_res_penalty_qual_obj','-append');
% 
%     % Load penalty cutoff 
% high_res_penalty_goals = load('test_save_files/high_resolution_penalty_goals.mat'); 
% penalty_cutoff = high_res_penalty_goals.penalty_goals; 
% 
% satisfy_flag_matrix_SSE = streamline_penalty_criteria(result_file_SSE.all_high_res_penalty_SSE,penalty_cutoff);
% satisfy_flag_matrix_qual_obj = streamline_penalty_criteria(result_file_qual_obj.all_high_res_penalty_qual_obj,penalty_cutoff);
% 
% num_satisfied_criteria_per_sample_streamlined_SSE = sum(satisfy_flag_matrix_SSE,2); 
% num_satisfied_criteria_per_sample_streamlined_qual_obj = sum(satisfy_flag_matrix_qual_obj,2); 
% num_satisfied_sample_all_SSE = nan(size(satisfy_flag_matrix_SSE,2),1);
% num_satisfied_sample_all_qual_obj = nan(size(satisfy_flag_matrix_qual_obj,2),1);
% for num_satisfied = 1:size(satisfy_flag_matrix_SSE,2)
% 
%     num_satisfied_sample_SSE = sum(num_satisfied_criteria_per_sample_streamlined_SSE >= num_satisfied); 
%     num_satisfied_sample_all_SSE(num_satisfied) = num_satisfied_sample_SSE; 
% 
%     num_satisfied_sample_qual_obj = sum(num_satisfied_criteria_per_sample_streamlined_qual_obj >= num_satisfied); 
%     num_satisfied_sample_all_qual_obj(num_satisfied) = num_satisfied_sample_qual_obj; 
% end
% % 
% % % labels = {'>= 5','>= 10','>= 15','>= 20','>= 25','>= 30','>= 32'};
% % % ,'Labels',labels
% % funnelchart(num_satisfied_sample_all_qual_obj)
% % funnelchart(num_satisfied_sample_all_SSE)
% 
% filter_idx_list_SSE = find(num_satisfied_criteria_per_sample_streamlined_SSE >= 8); 
% filter_idx_list_qual_obj = find(num_satisfied_criteria_per_sample_streamlined_qual_obj >= 8); 
% 
% filtered_high_res_penalty_SSE = result_file_SSE.all_high_res_penalty_SSE(filter_idx_list_SSE,:);
% filtered_high_res_penalty_qual_obj = result_file_qual_obj.all_high_res_penalty_qual_obj(filter_idx_list_qual_obj,:);
% 
% 
% 
% plot_qual_obj_dist(filtered_high_res_penalty_SSE,penalty_cutoff)
% plot_qual_obj_dist(filtered_high_res_penalty_qual_obj,penalty_cutoff)
% 
% 
% species_name_list = simFunction.Observables.Name;
% species_option_struct = struct('promotor_oi','sigma70_strong','conc_idx',1);
% species_option_struct.species_name_list = species_name_list; 
% 
% 
% for k = 1:length(filter_idx_list_qual_obj)
%     idx_oi = filter_idx_list_qual_obj(k);
%     estimated_params_oi = result_file_qual_obj.all_estimated_params(idx_oi,:);
%     [simulated_time,simulated_data] = simFunction(estimated_params_oi,tEnd,dosing_information,tStart:tEnd);
% 
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
%     % plot_simulated(simulated_time,simulated_data,'BundledSpecies',[],[],species_option_struct)
%     % plot_simulated(simulated_time,simulated_data,'exhaustive')
% end