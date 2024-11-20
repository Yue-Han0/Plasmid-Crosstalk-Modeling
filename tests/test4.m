clear
clc

currentpath = pwd; 
addpath(genpath(currentpath))
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))
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

%% Run missing in new mechanism sampling 
% load('test_save_files/202410_TestBindingSites_sampled.mat','all_sampled_params','all_high_res_penalty_terms','problemObject','dosing_information');
% num_missing = sum(all(isnan(all_high_res_penalty_terms),2)); 
% num_rerun = 0; 
% for iter = 1:size(all_high_res_penalty_terms,1)
%     high_res_penalty_term = all_high_res_penalty_terms(iter,:);
%     if all(isnan(high_res_penalty_term))
%         sampled_params = all_sampled_params(iter,:); 
%         [scaled_penalty_vec,~,~,~] = wrapper_calculate_obj_higher_resolution(sampled_params,problemObject,...
%             dosing_information,[],[],mode,num_conc,num_promotor);
%         all_high_res_penalty_terms(iter,:) = scaled_penalty_vec; 
%         num_rerun = num_rerun + 1; 
%     end
% 
% 
% end
for run_idx = 1:3841

    result_file_name = sprintf('param_est_run_save/20241007_TestBindingSites_param_sampling_constrained_run%d.mat',run_idx); 
    if ~exist(result_file_name,'file')
        run_param_and_init_cond_sampling_v2_NewMech(1,run_idx);
    end
end




%% Plot 1st-order k_deg vs. increasing mRNA concentration
% load('test_save_files/20240912_fitted_firstOrderRNAdeg_k.mat')
% init_mRNA_conc = output_table.initial_mRNA_conc;
% k_deg_1 = output_table.k_deg_1;
% k_deg_2 = output_table.k_deg_2; 
% k_deg_3 = output_table.k_deg_3;
% mean_k_deg = mean([k_deg_1 k_deg_2 k_deg_3],2);
% sd_k_deg = std([k_deg_1 k_deg_2 k_deg_3],0,2);
% 
% figure; 
% scatter(init_mRNA_conc(4:end),mean_k_deg(4:end),'filled'); 
% hold on
% for conc_idx = 4:length(init_mRNA_conc)
%     plot([init_mRNA_conc(conc_idx),init_mRNA_conc(conc_idx)],...
%         [mean_k_deg(conc_idx) - sd_k_deg(conc_idx),mean_k_deg(conc_idx) + sd_k_deg(conc_idx)],...
%         'LineWidth',1.5,'Color','k')
% 
% end
% set(gca,'FontSize',12)
% xlabel('Initial mRNA concentration(nM)')
% ylabel('Estimated k_{deg} (s^{-1})')

%% Compile penalty terms with high resolution 

% all_high_res_penalty_terms = nan(1591 * 100,32); 
% for run_idx = 1:1591
%     result_file_name = sprintf('param_est_run_save/20240821_param_sampling_constrained_run%d.mat',run_idx);
%     result_file = load(result_file_name);
% 
%     if isequal(run_idx,1)
%         problemObject = result_file.problemObject;
%         dosing_information = resultsave(sprintf('param_est_run_save/20240821_param_sampling_constrained_run%d.mat',run_idx),'high_res_all_penalty_terms','high_res_all_penalty_terms_unscaled',...
%     'sampled_params_selected','dosing_information');_file.dosing_information;
%     end
% 
%     high_res_all_penalty_terms = result_file.high_res_all_penalty_terms; 
%     all_high_res_penalty_terms((run_idx - 1) * 100 + 1:run_idx * 100,:) = high_res_all_penalty_terms;
% 
% end
% 
% save('test_save_files/202409_high_res_penalty_sampling_result_summary.mat','all_high_res_penalty_terms')

%% Run missing samples 
% missing_run_idx_list = []; 
% for run_idx = 1:1591
%     result_file_name = sprintf('param_est_run_save/20240821_param_sampling_constrained_run%d.mat',run_idx);
% 
%     if ~exist(result_file_name,'file')
%         missing_run_idx_list = [missing_run_idx_list,run_idx];
% 
%     end
% end
% 
% for j = 1:length(missing_run_idx_list)
%     missing_run_idx= missing_run_idx_list(j); 
%     fprintf('Running iter #%d',missing_run_idx)
%     run_param_and_init_cond_sampling_v2(missing_run_idx)
% end


%% Review TX-level crosstalk 
% TXcrosstalk_file_name = 'param_est_run_save/20230509_param_est_run1509_.mat'; 
% TXcrosstalk_file = load(TXcrosstalk_file_name); 
% 
% TX_problemObject_updated = TXcrosstalk_file.problemObject; 
% [~,min_idx] = min(TXcrosstalk_file.metric_summary); 
% opt_fitResult = TXcrosstalk_file.all_fitResults{min_idx}; 
% opt_estimated_params = [opt_fitResult.ParameterEstimates.Estimate];
% estimated_params_TX = {TX_problemObject_updated.Estimated.Name};
% 
% PE_sampling_result_file_name = '20240717_sampling_constrained_result_summary.mat';
% PE_samping_result_file = load(PE_sampling_result_file_name); 
% 
    % Check the value of estimated TX parameters in the context of PE
    % parameter sampling 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject_ctrl = setProblemObject_v2(group_description,param_info_path,[]); 
% estimated_params_PE = {problemObject_ctrl.Estimated.Name}; 
% 
% all_estimated_params_temp = [estimated_params_TX,estimated_params_PE];
% all_estimated_params = unique(all_estimated_params_temp);
% 
% figure; 
% for param_idx = 1:length(all_estimated_params)
%     subplot(ceil(sqrt(length(all_estimated_params))),ceil(sqrt(length(all_estimated_params))),param_idx)
% 
%     estimated_params_name = all_estimated_params{param_idx}; 
% 
%     % Find if estimated param in TX & its value 
%     estimated_param_idx_TX = strcmp(estimated_params_TX,estimated_params_name); 
%     if any(estimated_param_idx_TX)
%         estimated_param_val_TX = opt_estimated_params(estimated_param_idx_TX);
%     else
%         estimated_param_val_TX = nan; 
%     end
% 
%     % Find if estiamted param in PE & its value 
%     estimated_param_idx_PE = strcmp(estimated_params_PE,estimated_params_name); 
%     if any(estimated_param_idx_PE)
%         sampled_param_val_PE = PE_samping_result_file.all_sampled_params(:,estimated_param_idx_PE);
%     else
%         sampled_param_val_PE = nan; 
%     end
% 
%     % Plot both value
%     if all(~isnan(sampled_param_val_PE))
%         histogram(log10(sampled_param_val_PE));
%     end
%     hold on 
%     if ~isnan(estimated_param_val_TX)
%         xline(log10(estimated_param_val_TX),'LineWidth',1.5,'Color','r'); 
%     end
%     title(strrep(estimated_params_name,'_',' ')); 
% 
% end

%     % Fit TX-level crosstalk data using updated model 
% group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% param_info_path = 'parameters_v3.xlsx';
% TX_problemObject_updated = setProblemObject_v3(group_description,[],param_info_path); 
% save_path_id = 'TX_updated'; 
% wrapper_fit_to_data_w_probObject_v2(TX_problemObject_updated,save_path_id)
% 
%     % Check the amount of mRNA produced in fitted TX-level crosstalk model
% model_species_name_list = {TX_problemObject.Model.Species.Name};
% species_oi_name = 'RNA utrbroc'; 
% species_oi_idx = strcmp(model_species_name_list,species_oi_name); 
% additional_track_species = [model_species_name_list(1:species_oi_idx - 1),model_species_name_list(species_oi_idx + 1:end)];
% simFunction_TX = create_simFun_from_problemObject(TX_problemObject,additional_track_species);
% dosing_information = create_dosing_info_from_problemObject(TX_problemObject); 
% 
% [simulated_time,simulated_data] = simFunction(opt_estimated_params,tEnd,dosing_information,tStart:tEnd);
% 
% 
    % Maybe use the sampled parameters in PE space and check phenotypes in
    % TX-level crosstalk? 

%% Troubleshoot weights in fitProblem struct 
% group_description = {'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject = setProblemObject_v3_normalized(group_description,[],param_info_path);
% problemObject_ctrl = setProblemObject_v3(group_description,[],param_info_path);
% save_path_id = 'weight_on_subset'; 

    % Confirm that weights are actually applied correctly
% 
% result_fileName_1_weights = 'param_est_run_save/20230812_param_est_run0756_weight_on_subset.mat';
% result_fileName_2_no_weights = 'param_est_run_save/20230812_param_est_run1942_weight_on_subset.mat';
% result_file_1_weights = load(result_fileName_1_weights);
% result_file_2_no_weights = load(result_fileName_2_no_weights);
% 
% wrapper_analyze_fit_result_v2(result_file_1.all_fitResults,result_file_1.problemObject,...
%     result_fileName_1,'SSE',false,false); 
% wrapper_analyze_fit_result_v2(result_file_2.all_fitResults,result_file_2.problemObject,...
%     result_fileName_2,'SSE',false,false); 
% 
% 
%     % Check if problemObject contains the same species, reaction, and
%     % parameters 
% result_fileName_weights = 'param_est_run_save/20230719_param_est_run1343_weight_on_subset.mat'; 
% result_file_weights = load(result_fileName_weights);
% problemObject_ori_weights = result_file_weights.problemObject; 
% 
% wrapper_fit_to_data_w_probObject(problemObject,save_path_id)
% 
% Directly compare fitResult with and without added weights 
% result_fileName_weights = 'param_est_run_save/20230809_param_est_run1456_weight_on_subset.mat'; 
% result_fileNmae_ctrl_no_weights = 'param_est_run_save/20230809_param_est_run1105_no_weight_on_subset_ctrl.mat'; 
% 
% result_file_weights = load(result_fileName_weights); 
% result_file_ctrl_no_weights = load(result_fileNmae_ctrl_no_weights);
% 
% Compare optFitResult struct 
% [~,sort_idx_weights] = sort(result_file_weights.metric_summary,'ascend');
% [~,sort_idx_ctrl_no_weights] = sort(result_file_ctrl_no_weights.metric_summary,'ascend'); 
% opt_fitResult_weights = result_file_weights.all_fitResults{sort_idx_weights(1)};
% opt_fitResult_ctrl_no_weights = result_file_ctrl_no_weights.all_fitResults{sort_idx_ctrl_no_weights(1)}; 
% 
% % Try to verify the reported metric 
% ctrl_no_weight_verify_SSE = sum(opt_fitResult_ctrl_no_weights.R.^2); 
% weight_verify_SSE = sum(opt_fitResult_weights.R .^2); 
% 
% wrapper_analyze_fit_result_v2(result_file_ctrl_no_weights.all_fitResults,result_file_ctrl_no_weights.problemObject,...
%     result_fileNmae_ctrl_no_weights,'SSE',false,false); 
% 
% figure;
% histogram(log10(result_file_1_weights.metric_summary))
% hold on 
% histogram(log10(result_file_2_no_weights.metric_summary))
% legend('w\ weights','wo\ weights')
% 
% % What about capturing the phenotypes? 
%     % Calculating qualitative penalties and comparing distributions 
% all_crosstalk_ratios_weight = nan(length(result_file_1_weights.all_fitResults),num_conc);
% all_crosstalk_ratios_no_weight = nan(length(result_file_2_no_weights.all_fitResults),num_conc);
% all_estimated_params_weight = nan(length(result_file_1_weights.all_fitResults),length(result_file_1_weights.problemObject.Estimated)); 
% all_estimated_params_no_weight = nan(length(result_file_2_no_weights.all_fitResults),length(result_file_2_no_weights.problemObject.Estimated)); 
% 
% 
%     % Extract problemObject and create simFunction
% problemObject_weight = result_file_1_weights.problemObject;
% problemObject_no_weight = result_file_2_no_weights.problemObject;
%         % Reorganize observables such that sfGFP comes first 
% model = problemObject_weight.Model; 
% species_name_list = {model.Species.Name}; 
% remove_idx = strcmp('protein sfGFP*',species_name_list); 
% additional_track_species_list = species_name_list(~remove_idx); 
% observables = [{'protein sfGFP*'},additional_track_species_list];
% simFunction_weight = create_simFun_from_problemObject(problemObject_weight,additional_track_species_list); 
% simFunction_no_weight = create_simFun_from_problemObject(problemObject_no_weight,additional_track_species_list); 
% 
% dosing_information = create_dosing_info_from_problemObject(problemObject_weight);
% sampled_param_names = {problemObject_weight.Estimated.Name};
% target_name_list_struct = load('test_save_files/20240717_param_init_cond_sampling_constrained.mat','target_name_list');
% target_name_list = target_name_list_struct.target_name_list; 
% 
% for iter = 1:length(result_file_1_weights.all_fitResults)
% 
%     % Get estimated params 
%     fitResult_oi_weight = result_file_1_weights.all_fitResults{iter}; 
%     fitResult_oi_no_weight = result_file_2_no_weights.all_fitResults{iter}; 
%     estimated_params_weight = fitResult_oi_weight.ParameterEstimates.Estimate; 
%     estimated_params_no_weight = fitResult_oi_no_weight.ParameterEstimates.Estimate; 
% 
%     % Calculate crosstalk ratios for each case 
%     [simulated_time_weight,simulated_data_weight] = simFunction_weight(estimated_params_weight',tEnd,dosing_information,tStart:tEnd);
%     [simulated_time_no_weight,simulated_data_no_weight] = simFunction_no_weight(estimated_params_no_weight',tEnd,dosing_information,tStart:tEnd);
%     crosstalk_ratios_weight = nan(1,num_conc);
%     crosstalk_ratios_no_weight = nan(1,num_conc);
%     for conc_idx = 1:num_conc
%         crosstalk_ratios_weight(conc_idx) = max(simulated_data_weight{conc_idx + num_conc}(:,1)) ./ max(simulated_data_weight{conc_idx}(:,1)); 
%         crosstalk_ratios_no_weight(conc_idx) = max(simulated_data_no_weight{conc_idx + num_conc}(:,1)) ./ max(simulated_data_no_weight{conc_idx}(:,1)); 
%     end
% 
%     % Allocate 
%     all_crosstalk_ratios_weight(iter,:) = crosstalk_ratios_weight; 
%     all_crosstalk_ratios_no_weight(iter,:) = crosstalk_ratios_no_weight; 
%     all_estimated_params_weight(iter,:) = estimated_params_weight; 
%     all_estimated_params_no_weight(iter,:) = estimated_params_no_weight;
% 
% end
% 
% % Plot out crosstalk ratio for all 48? 
% figure; 
% for iter = 1:num_iter
%     subplot(ceil(sqrt(num_iter)),ceil(sqrt(num_iter)),iter)
% 
%     plot(conc_vec,all_crosstalk_ratios_weight(iter,:),'LineWidth',1.5,'Color','r')
%     hold on
%     plot(conc_vec,all_crosstalk_ratios_no_weight(iter,:),'LineWidth',1.5,'Color','k')
%     title(sprintf('Iter #%d',iter))
% 
% 
% end
% legend('w\ weights','wo\ weights')
% sgtitle('Crosstalk Ratio Comparison - weights in data fitting')
% 
% % Count # cases where negative/positive crosstalk/both are captured in each
% num_neg_crosstalk_captured = zeros(1,2); % 1st column - with weights; 2nd column - without weights 
% num_pos_crosstalk_captured = zeros(1,2); 
% num_transition_captured = zeros(1,2);
% for iter = 1:num_iter
%     neg_crosstalk_flag_1 = any(all_crosstalk_ratios_weight(iter,4:7) < 0.9); 
%     neg_crosstalk_flag_2 = any(all_crosstalk_ratios_no_weight(iter,4:7) < 0.9); 
% 
%     pos_crosstalk_flag_1 = any(all_crosstalk_ratios_weight(iter,1:4) > 1.1); 
%     pos_crosstalk_flag_2 = any(all_crosstalk_ratios_no_weight(iter,1:4) > 1.1); 
% 
%     num_neg_crosstalk_captured(1) = num_neg_crosstalk_captured(1) + neg_crosstalk_flag_1;
%     num_neg_crosstalk_captured(2) = num_neg_crosstalk_captured(2) + neg_crosstalk_flag_2;
%     num_pos_crosstalk_captured(1) = num_pos_crosstalk_captured(1) + pos_crosstalk_flag_1; 
%     num_pos_crosstalk_captured(2) = num_pos_crosstalk_captured(2) + pos_crosstalk_flag_2; 
% 
%     num_transition_captured(1) = num_transition_captured(1) + (neg_crosstalk_flag_1 && pos_crosstalk_flag_1); 
%     num_transition_captured(2) = num_transition_captured(2) + (neg_crosstalk_flag_2 && pos_crosstalk_flag_2); 
% 
% end
% 

% result_fileName_3highConc = 'param_est_run_save/20230717_param_est_run1557_sigma70_weak_3_high_conc.mat'; 

% result_fileName_weights = 'param_est_run_save/20230719_param_est_run1343_weight_on_subset.mat'; 

% result_fileName_weights = 'param_est_run_save/202'
%% Check parameter assignment in RNA degradation 

% keys = {'TX','PE'};
% values = {[],[89:91,96:98]}; 
% group_num = dictionary(keys,values);
% 
% param_info_path = 'parameters_v3_1.xlsx';
% all_save_path_id_list = {'Ctrl_ori_params','Ctrl_exp_params','PolysomeOnly_ori_params','PolySomeOnly_exp_params','NonProcDegOnly_ori_params',...
%     'NonProcDegOnly_exp_params','PolySomeNNonProcDeg_ori_params','PolySomeNNonProcDeg_exp_params'}; 
% 
% init_problemObject = setProblemObject_v3([],group_num,param_info_path);
% 
% % problemObject_PolySomeOnly = modify_mechanism_in_problemObject(init_problemObject,'PolySomeOnly'); 
% % problemObject_NonProcDegOnly = modify_mechanism_in_problemObject(init_problemObject,'NonProcDegOnly',param_info_path); 
% problemObject_PolySomeNNonProcDeg = modify_mechanism_in_problemObject(init_problemObject,'PolySomeNNonProcDeg',param_info_path); 
% 
% problemObject = problemObject_PolySomeNNonProcDeg; 
% % For each estimated parameters in RNA degradation, find the reaction that
% % uses the parameter and verify validity 
% estimated_param_name_list = {problemObject.Estimated.Name};
% model_rxns = problemObject.Model.Reactions; 
% model_rxn_rate = {problemObject.Model.Reactions.ReactionRate}; 
%     % Construct a dictionary whose key is parameter name and value is
%     % reaction 
% keys_param_names = {};
% values_reaction_names = {}; 
% values_reaction_rates = {}; 
% for est_idx = 1:length(estimated_param_name_list)
%     estimated_param_name = estimated_param_name_list{est_idx};
%     if contains(estimated_param_name,'RNAdeg') % Check for each RNA deg related parameter
%         keys_param_names{end + 1} = estimated_param_name; 
%         relevant_rxn_idx_list = contains(model_rxn_rate,strcat(estimated_param_name ,'*')) | contains(model_rxn_rate,strcat('*',estimated_param_name ))...
%             | contains(model_rxn_rate,estimated_param_name  + " ") | contains(model_rxn_rate," " + estimated_param_name ); 
%         values_reaction_names{end + 1} = model_rxns(relevant_rxn_idx_list); 
%         values_reaction_rates{end + 1} = model_rxn_rate(relevant_rxn_idx_list); 
%     end
% end
% RNAdegParam_in_Rxn_dict = dictionary(keys_param_names,values_reaction_names); 



%% Analyze fitting results with polysome and new RNA degradation mechanism 
% param_option_list = {'ori','exp'};
% mechanism_option_list = {'Ctrl','PolysomeOnly','NonProcDegOnly','PolySomeNNonProcDeg'};
% 
% %     % Load result files 
% %         % Get all filenames in the directory 
% % directory_path = 'param_est_run_save'; 
% % directory_files = dir(directory_path); 
% % file_names = {directory_files.name};
% % for param_option_idx = 1:length(param_option_list)
% %     param_option = param_option_list{param_option_idx};
% %     for mech_idx = 1:length(mechanism_option_list)
% %         mechanism_option = mechanism_option_list{mech_idx}; 
% %         variable_name_suffix = sprintf('_%s_%s',mechanism_option,param_option);
% %         result_fileName_idx = contains(file_names,param_option) & contains(file_names,mechanism_option); 
% %         result_fileName = file_names{result_fileName_idx}; 
% %         result_file = load(result_fileName);
% %         % wrapper_analyze_fit_result_v2(result_file.all_fitResults,result_file.problemObject,result_fileName,'SSE',false,false); 
% %         eval(sprintf('all_fitResults_%s = result_file.all_fitResults;',variable_name_suffix)); 
% %         eval(sprintf('problemObject_%s = result_file.problemObject;',variable_name_suffix)); 
% %         eval(sprintf('metric_summary_%s = result_file.metric_summary;',variable_name_suffix)); 
% %     end
% % end
% 
%     % Visualize SSE in each case 
% figure; 
% for mech_idx = 1:length(mechanism_option_list)
%     mechanism_option = mechanism_option_list{mech_idx}; 
%     for param_option_idx = 1:length(param_option_list)
%         param_option = param_option_list{param_option_idx};
%         variable_name_suffix = sprintf('_%s_%s',mechanism_option,param_option);
% 
%         subplot(4,2,(mech_idx - 1) * length(param_option_list) + param_option_idx)
% 
%         % Get SSE distribution and plot in historgram 
%         eval(sprintf('metric_summary_oi = metric_summary_%s;',variable_name_suffix)); 
%         eval(sprintf('[~,sorted_SSE_idx_%s] = sort(metric_summary_%s,''ascend'');',...
%             variable_name_suffix,variable_name_suffix)); 
%         histogram(log10(metric_summary_oi)); 
%         xlim([6,9])
% 
%         % Report number of searched out of 48 where large positive
%         % crosstalk is captured 
%             % Check positive crosstalk ratio for the lowest reporter
%             % plasmid concentration 
%         num_large_positive_crosstalk = 0; 
%         num_positive_crosstalk = 0; 
%         for iter = 1:num_iter
%             eval(sprintf('fitResult_oi = all_fitResults_%s{iter};',variable_name_suffix));
%             simData_oi = fitted(fitResult_oi); 
%             simData_opt_single_no_empty_lowest = simData_oi(1); 
%             simData_opt_single_empty_lowest = simData_oi(length(simData_oi)/2 + 1); 
%             species_oi_idx = strcmp(simData_opt_single_no_empty_lowest.DataNames,species_oi);
%             approx_positive_crosstalk_ratio = max(simData_opt_single_empty_lowest.Data(:,species_oi_idx)) ./ ...
%                 max(simData_opt_single_no_empty_lowest.Data(:,species_oi_idx)); 
% 
%             if approx_positive_crosstalk_ratio > 1.1 % only plot time course if capturing positive crosstalk 
%                 num_positive_crosstalk = num_positive_crosstalk + 1; 
%             end
%             if approx_positive_crosstalk_ratio > 2
%                 num_large_positive_crosstalk = num_large_positive_crosstalk + 1; 
%             end
%         end
%         title(sprintf('%s \n %d out of 48 captured positive crosstalk \n %d out of 48 captured large positive crosstalk',...
%             strrep(variable_name_suffix,'_',' '),num_positive_crosstalk,num_large_positive_crosstalk)); 
%     end
% end
% 
%     % For each scenario, plot out no empty vs. empty & exp vs. fitted 
%     % For each result file, select top 5 best fitted model and check crosstalk
%     % ratio 
% 
% % species_oi = 'protein sfGFP*'; 
% % 
% % % Get experimental data 
% % experimental_groupedData = problemObject__Ctrl_ori.Data; 
% % 
% % for mech_idx = 1:length(mechanism_option_list)
% %     mechanism_option = mechanism_option_list{mech_idx}; 
% %     for param_option_idx = 1:length(param_option_list)
% %         param_option = param_option_list{param_option_idx};
% %         variable_name_suffix = sprintf('_%s_%s',mechanism_option,param_option);
% % 
% %         figure; 
% %         for idx = 1:5
% %             % Select fitResult of interest and get simData
% %             eval(sprintf('fitResult_oi = all_fitResults_%s{sorted_SSE_idx_%s(idx)};',variable_name_suffix,variable_name_suffix)); 
% %             simData_oi = fitted(fitResult_oi);
% %             num_conc = length(simData_oi)/2; 
% % 
% %             for conc_idx = 1:num_conc
% % 
% %                 subplot(5,num_conc,(idx - 1) * num_conc + conc_idx)
% % 
% %                 % Plot simulated data 
% %                 simData_no_empty_oi = simData_oi(conc_idx); 
% %                 simData_empty_oi = simData_oi(num_conc + conc_idx); 
% %                 species_oi_idx = find(strcmp(simData_no_empty_oi.DataNames,species_oi)); 
% %                 plot(simData_no_empty_oi.Time,simData_no_empty_oi.Data(:,species_oi_idx),'--','LineWidth',1.5,'Color','g')
% %                 hold on 
% %                 plot(simData_empty_oi.Time,simData_empty_oi.Data(:,species_oi_idx),'--','LineWidth',1.5,'Color',[0.5,0.5,0.5])
% % 
% %                 % Plot experimental data 
% %                 group_num_no_empty = 168 + conc_idx; 
% %                 group_num_empty = 168 + 7 + conc_idx; 
% %                 no_empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_no_empty); 
% %                 no_empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_no_empty); 
% %                 empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_empty); 
% %                 empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_empty); 
% %                 plot(no_empty_timeVec,no_empty_GFP_conc,'Color','g','LineWidth',1.5)
% %                 hold on
% %                 plot(empty_timeVec,empty_GFP_conc,'Color',[0.5,0.5,0.5],'LineWidth',1.5)
% % 
% %                 title(sprintf('%.1f nM Plasmid',conc_vec(conc_idx))); 
% %             end
% % 
% %             if isequal(conc_idx,1)
% %                 ylabel(sprintf('#%d Best Fitting',idx))
% %             end
% %         end
% %         legend('Sim - no empty','Sim - empty','Exp - no empty','Exp - empty')
% %         sgtitle(strrep(variable_name_suffix,'_',' '))
% % 
% %     end
% % end
% 
% %     % For cases where positive crosstalk is captured, check the residue
% %     % mRNA concentration (& other species concentrations) 
% %         % For control and Polysome Only, only check the best fitting 
% % for mech_idx = 1:2%length(mechanism_option_list)
% %     mechanism_option = mechanism_option_list{mech_idx}; 
% %     for param_option_idx = 1:length(param_option_list)
% %         param_option = param_option_list{param_option_idx};
% %         variable_name_suffix = sprintf('_%s_%s',mechanism_option,param_option);
% %         eval(sprintf('fitResult_opt = all_fitResults_%s{sorted_SSE_idx_%s(1)}',variable_name_suffix,variable_name_suffix)); 
% %         simData_opt = fitted(fitResult_opt);
% %         num_conc = length(simData_opt)/2; 
% %             % Let's just plot out resource utilization by the largest
% %             % reporter plasmid concentration fitted 
% %         simData_opt_single_no_empty = simData_opt(num_conc); 
% %         simData_opt_single_empty = simData_opt(end); 
% %         all_species_name_list = simData_opt_single_no_empty.DataNames; 
% %         figure; 
% %         for species_idx = 1:length(all_species_name_list)
% %             subplot(ceil(sqrt(length(all_species_name_list))),ceil(sqrt(length(all_species_name_list))),species_idx)
% %             plot(simData_opt_single_no_empty.Time,simData_opt_single_no_empty.Data(:,species_idx),'LineWidth',1.5,'Color','g')
% %             hold on 
% %             plot(simData_opt_single_empty.Time,simData_opt_single_empty.Data(:,species_idx),'LineWidth',1.5,'Color',[0.5,0.5,0.5])
% %             species_name_for_plot = strrep(all_species_name_list{species_idx},'TXTL','');
% %             species_name_for_plot = strrep(species_name_for_plot,'_',' ');
% %             title(species_name_for_plot)
% %         end
% %         sgtitle(strrep(variable_name_suffix,'_',' ')); 
% %     end
% % end
% 
% %     % For Non-processive RNA deg, check each case positive crosstalk is
% %     % captured 
% % species_oi_list = {'RNase','RNA utrGFP--sfGFP','RNA utrkanR--kanR','RNase_bound','mRNA_deact'}; 
% % for mech_idx = 3:length(mechanism_option_list)
% %     mechanism_option = mechanism_option_list{mech_idx}; 
% %     for param_option_idx = 1:length(param_option_list)
% %         param_option = param_option_list{param_option_idx};
% %         variable_name_suffix = sprintf('_%s_%s',mechanism_option,param_option);
% % 
% %         figure; 
% %         for iter = 1:10
% % 
% %             eval(sprintf('fitResult_oi = all_fitResults_%s{sorted_SSE_idx_%s(iter)};',variable_name_suffix,variable_name_suffix)); 
% %             simData_opt = fitted(fitResult_oi);
% %             num_conc = length(simData_opt)/2; 
% % 
% %                 % Check positive crosstalk ratio for the lowest reporter
% %                 % plasmid concentration 
% %             simData_opt_single_no_empty_lowest = simData_opt(1); 
% %             simData_opt_single_empty_lowest = simData_opt(num_conc + 1); 
% %             species_oi_idx = strcmp(simData_opt_single_no_empty_lowest.DataNames,species_oi);
% %             approx_positive_crosstalk_ratio = max(simData_opt_single_empty_lowest.Data(:,species_oi_idx)) ./ ...
% %                 max(simData_opt_single_no_empty_lowest.Data(:,species_oi_idx)); 
% % 
% %             if approx_positive_crosstalk_ratio > 1.1 % only plot time course if capturing positive crosstalk 
% %                     % Plot out time course for species of interest for the
% %                     % highest reporter plasmid concentrations fitted 
% %                 simData_opt_single_no_empty = simData_opt(num_conc); 
% %                 simData_opt_single_empty = simData_opt(end); 
% %                 for species_idx = 1:length(species_oi_list)
% %                     species_oi_name = species_oi_list{species_idx}; 
% %                     species_oi_idx = find(strcmp(simData_opt_single_no_empty.DataNames,species_oi_name)); 
% %                     subplot(10,length(species_oi_list),(iter - 1) * length(species_oi_list) + species_idx)
% %                     plot(simData_opt_single_no_empty.Time,simData_opt_single_no_empty.Data(:,species_oi_idx),'LineWidth',1.5,'Color','g')
% %                     hold on 
% %                     plot(simData_opt_single_empty.Time,simData_opt_single_empty.Data(:,species_oi_idx),'LineWidth',1.5,'Color',[0.5,0.5,0.5])
% %                     species_name_for_plot = strrep(species_oi_name,'_',' ');
% %                     title(species_name_for_plot)
% %                 end
% % 
% %             end
% %         end
% %         sgtitle(strrep(variable_name_suffix,'_',' ')); 
% %     end
% % end
% 
% % Manually modify parameter values for the non-processive RNA degradation
% % mechanism to explore if RNA can be fully degraded without compromising
% % captured positive crosstalk
% 
%     % For non-processive with original parameter range, compile estimated
%     % parameters 
% 
% num_params_NonProcDegOnly = length(problemObject__NonProcDegOnly_ori.Estimated); 
% all_estimated_params = nan(num_iter,num_params_NonProcDegOnly); 
% pos_crosstalk_flag = nan(num_iter,1); 
% species_oi = 'protein sfGFP*';
% % 
% % for iter = 1:num_iter
% %     fitResult = all_fitResults__NonProcDegOnly_ori{iter}; 
% %     simData = fitted(fitResult); 
% % 
% %     % Check positive crosstalk ratio for the lowest reporter
% %     % plasmid concentration 
% %     simData_single_no_empty_lowest = simData(1); 
% %     simData_single_empty_lowest = simData(4); 
% %     species_oi_idx = strcmp(simData_single_no_empty_lowest.DataNames,species_oi);
% %     approx_positive_crosstalk_ratio = max(simData_single_empty_lowest.Data(:,species_oi_idx)) ./ ...
% %         max(simData_single_no_empty_lowest.Data(:,species_oi_idx));     
% %     pos_crosstalk_flag(iter) = approx_positive_crosstalk_ratio > 1.1; 
% % 
% %     % Compile parameters 
% %     all_estimated_params(iter,:) = [fitResult.ParameterEstimates.Estimate]; 
% % end
% % 
% %     % Check parameter distribution for cases capturing positive crosstalk
% %     % and cases that do not 
% % pos_crosstalk_captured_idx_list = find(pos_crosstalk_flag); 
% % pos_crosstalk_not_captured_idx_list = find(~pos_crosstalk_flag); 
% % estimated_params_captured = all_estimated_params(pos_crosstalk_captured_idx_list,:); 
% % estimated_params_not_captured = all_estimated_params(pos_crosstalk_not_captured_idx_list,:);
% % 
% % estimated_param_names = {problemObject__NonProcDegOnly_ori.Estimated.Name}; 
% % for param_idx = 1:length(estimated_param_names)
% %     param_name = estimated_param_names{param_idx}; 
% %     subplot(ceil(sqrt(length(estimated_param_names))),ceil(sqrt(length(estimated_param_names))),param_idx)
% %     histogram(log(estimated_params_captured(:,param_idx)))
% %     hold on 
% %     histogram(log(estimated_params_not_captured(:,param_idx)))
% %     title(strrep(param_name,'_',' '));
% % 
% % end
% 
%     % Select the best fitted case and modify parameter values 
%         % Get the best fitting parameters 
% [~,sort_idx_NonProcDegOnly] = sort(metric_summary__NonProcDegOnly,'ascend'); 
% fitResult_opt = all_fitResults__NonProcDegOnly{sort_idx_NonProcDegOnly(1)}; 
% estimated_params = [fitResult_opt.ParameterEstimates.Estimate]; 
% estimated_param_names = {fitResult_opt.ParameterEstimates.Name};
%         % Create simFunction 
% all_species_name = {problemObject__NonProcDegOnly_ori.Model.Species.Name};
% species_oi_idx = find(strcmp(all_species_name,species_oi)); 
% additional_tracked_species = [{'protein sfGFP*'},all_species_name(1:species_oi_idx-1),all_species_name(species_oi_idx+1:end)]; 
% simFunction = create_simFun_from_problemObject(problemObject__NonProcDegOnly_ori); 
% dosing_information = create_dosing_info_from_problemObject(problemObject__NonProcDegOnly_ori); 
%         % Modify parameters 
% modified_params = estimated_params;
% modify_param_name_list = {'TXTL_UTR_GFP_R','TXTL_UTR_kanR_R','TXTL_TL_init_k'}; 
% modify_target_value_list = [1e-09,1e-04,1e+04]; 
% for modify_idx = 1:length(modify_param_name_list)
%     modify_param_name = modify_param_name_list{modify_idx}; 
%     param_idx = find(strcmp(modify_param_name,estimated_param_names));
%     modified_params(param_idx) = modify_target_value_list(modify_idx); 
% end
% [simulated_time,simulated_data] = simFunction(modified_params,tEnd,dosing_information,tStart:tEnd);

%% Update parameters_v3.xlsx
% estimated_params_file = readtable('param_info/parameters_v3.xlsx'); 
% estimated_params_names = estimated_params_file.Name;
% 
% % Fix RNA degradation parameters 
% mRNA_params_file = load('test_save_files/20240530_mRNA_degradation_param_bounds.mat'); 
% mRNA_deg_params_CI = mRNA_params_file.param_bounds_CI; 
% mRNA_deg_params_names = mRNA_params_file.estimated_params_labels; 

% for mRNA_deg_param_idx = 1:length(mRNA_deg_params_names)
%     mRNA_deg_param_name = mRNA_deg_params_names{mRNA_deg_param_idx}; 
%     estimated_idx = strcmp(mRNA_deg_param_name,estimated_params_names); 
%     estimated_params_file{estimated_idx,'LB'} = mRNA_deg_params_CI(mRNA_deg_param_idx,1); 
%     estimated_params_file{estimated_idx,'UB'} = mRNA_deg_params_CI(mRNA_deg_param_idx,2); 
% 
%     % Update initial value
%     estimated_params_file{estimated_idx,'InitVal'} = mRNA_deg_params_CI(mRNA_deg_param_idx,1) + rand * (mRNA_deg_params_CI(mRNA_deg_param_idx,2) - mRNA_deg_params_CI(mRNA_deg_param_idx,1)); 
% 
% end
% 
% write(estimated_params_file,'param_info/parameters_v3.xlsx'); 

%% Calculate penalty terms at a high resolution
% result_file = load('test_save_files/20240701_sampling_constrained_result_summary.mat'); 
% all_sampled_params = result_file.all_sampled_params; 
% all_penalty_terms = result_file.all_calculated_penalty_terms; 
% sampled_params_file = load('test_save_files/20240613_param_init_cond_sampling_constrained.mat'); 
% all_dosing_information = sampled_params_file.all_dosing_information; 
% problemObject = sampled_params_file.problemObject; 
% kinetic_param_names = sampled_params_file.kinetic_param_names;
% target_name_list = sampled_params_file.target_name_list; 
% 
% for sample_idx = 1:size(all_penalty_terms,1)
% 
%     sampled_params = all_sampled_params(sample_idx,:);
%     if mod(sample_idx,100) == 0 
%         dosing_information = all_dosing_information{100}; 
%     else
%         dosing_information = all_dosing_information{mod(sample_idx,100)}; 
%     end
% 
%     temp_problemObject = problemObject; 
% 
%     [scaled_penalty_vec,penalty_term_labels,penalty_term_length,unscaled_penalty_vec] = ...
%         wrapper_calculate_obj_higher_resolution(sampled_params,temp_problemObject,dosing_information,kinetic_param_names,target_name_list,mode,num_conc,num_promotor);
%     if ~exist('high_res_penalty_terms','var') % Preassign 
%         high_res_penalty_terms = nan(size(all_penalty_terms,1),length(scaled_penalty_vec));
%         high_res_penalty_terms_unscaled = nan(size(all_penalty_terms,1),length(scaled_penalty_vec));
%     end
%     high_res_penalty_terms(sample_idx,:) = scaled_penalty_vec; 
%     high_res_penalty_terms_unscaled(sample_idx,:) = unscaled_penalty_vec; 
% 
% end


%% Compare calculated penalty terms for duplicate vs. no duplicate reactions 
% all_penalty_term_w_dup = nan(1000 * 10,6); 
% all_penalty_term_wo_dup = nan(1000 * 10,6); 
% 
% for iter = 1:10
%     % Load both result file using model with and without duplicate
%     % reactions 
%     old_result_file_name = sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d.mat',iter); 
%     new_result_file_name = sprintf('param_est_run_save/20240717_param_and_init_cond_sampling_constrained_run%d.mat',iter); 
%     old_result_file = load(old_result_file_name);
%     new_result_file = load(new_result_file_name); 
% 
%     % Check whether sampled parameters are the same (sanity check) 
%     old_result_file_sampled_params = old_result_file.sampled_params_selected;
%     new_result_file_sampled_params = new_result_file.sampled_params_selected; 
%     same_params_flag = all(all(old_result_file_sampled_params == new_result_file_sampled_params)); 
% 
%     % Record calculated penalty terms and compare 
%     if same_params_flag
%         all_penalty_term_w_dup((iter - 1) * 1000 + 1:iter * 1000,:) = old_result_file.all_penalty_terms; 
%         all_penalty_term_wo_dup((iter - 1) * 1000 + 1:iter * 1000,:) = new_result_file.all_penalty_terms; 
%     end
% 
%  end
% 
%  % Visualize changes in penalty terms 
% for penalty_idx = 1:length(old_result_file.penalty_labels)
%     subplot(2,3,penalty_idx)
%     filtered_penalty_term_w_dup_single = rmoutliers(all_penalty_term_w_dup(:,penalty_idx)); 
%     histogram(filtered_penalty_term_w_dup_single)
%     hold on 
%     filtered_penalty_term_wo_dup_single = rmoutliers(all_penalty_term_wo_dup(:,penalty_idx)); 
%     histogram(filtered_penalty_term_wo_dup_single)
%     title(strrep(old_result_file.penalty_labels{penalty_idx},'_',' '))
% 
% 
% end


%% Fit combinations of lower and higher plasmid concentrations & Play with weights
% param_info_path = 'parameters_v2.xlsx';
% 
% %     % 4 low concentrations 
% % keys = {'TX','PE'};
% % values = {[],[85:88,92:95]}; 
% % group_number = dictionary(keys,values);
% % problemObject = setProblemObject_v3([],group_number,param_info_path);
% % save_path_id = 'sigma70_weak_4_low_conc'; 
% result_fileName_4lowConc = 'param_est_run_save/20230717_param_est_run1554_sigma70_weak_4_low_conc.mat'; 
% % result_file_4lowConc = load(result_fileName_4lowConc); 
% 
% %     % 3 high concentrations 
% % keys = {'TX','PE'};
% % values = {[],[89:91,96:98]}; 
% % group_number_high_conc = dictionary(keys,values); 
% % problemObject = setProblemObject_v3([],group_number_high_conc,param_info_path);
% % save_path_id = 'sigma70_weak_3_high_conc'; 
% result_fileName_3highConc = 'param_est_run_save/20230717_param_est_run1557_sigma70_weak_3_high_conc.mat'; 
% % result_file_3highConc = load(result_fileName_3highConc);
% 
% %     % Weights on data 
% % group_description = {'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty'};
% % problemObject = setProblemObject_v3_normalized(group_description,[],param_info_path);
% % save_path_id = 'weight_on_subset'; 
% % 
% % wrapper_fit_to_data_w_probObject(problemObject,save_path_id)
% result_fileName_weights = 'param_est_run_save/20230719_param_est_run1343_weight_on_subset.mat'; 
% % result_file_weights = load(result_fileName_weights); 
% 
% %     % Analyze fitted results 
% % wrapper_analyze_fit_result_v2(result_file_4lowConc.all_fitResults,result_file_4lowConc.problemObject,result_fileName_4lowConc,'SSE',true,true)
% % wrapper_analyze_fit_result_v2(result_file_3highConc.all_fitResults,result_file_3highConc.problemObject,result_fileName_3highConc,'SSE',true,true)
% % wrapper_analyze_fit_result_v2(result_file_weights.all_fitResults,result_file_weights.problemObject,result_fileName_weights,'SSE',true,true)
% % 
% % % Plot metric distribution in each case
% %     % 4 low conc 
% % figure; 
% % histogram(log10(result_file_4lowConc.metric_summary)) 
% % xlabel('log10(SSE)')
% % ylabel('# Occurence')
% % title('Metric Distribution - fitting sigma70 weak 4 low conc')
% % 
% %     % 3 high conc
% % figure;
% % histogram(log10(result_file_3highConc.metric_summary)) 
% % xlabel('log10(SSE)')
% % ylabel('# Occurence')
% % title('Metric Distribution - fitting sigma70 weak 3 high conc')
% % 
% %     % weights 
% % figure;
% % histogram(log10(result_file_weights.metric_summary)) 
% % xlabel('log10(SSE)')
% % ylabel('# Occurence')
% % title('Metric Distribution - fitting sigma70 weak all conc scaled data')
% 
% species_name_oi = 'protein sfGFP*'; 
% experimental_groupedData = result_file_weights.problemObject.Data; 
% 
% % Plot simulated crosstalk in each case 
% [~,sort_idx_4lowConc] = sort(result_file_4lowConc.metric_summary,'ascend'); 
% [~,sort_idx_3highConc] = sort(result_file_3highConc.metric_summary,'ascend'); 
% [~,sort_idx_weights] = sort(result_file_weights.metric_summary,'ascend'); 
% 
% for best_idx = 1:10
% 
%     % Select fitResult of interest 
%     fitResult_oi_4lowConc = result_file_4lowConc.all_fitResults{sort_idx_4lowConc(best_idx)};  
%     fitResult_oi_3highConc = result_file_3highConc.all_fitResults{sort_idx_3highConc(best_idx)}; 
%     fitResult_oi_weights = result_file_weights.all_fitResults{sort_idx_weights(best_idx)}; 
% 
%     % Get fitted data 
%     fitted_data_4lowConc = fitted(fitResult_oi_4lowConc);
%     fitted_data_3highConc = fitted(fitResult_oi_3highConc);
%     fitted_data_weights = fitted(fitResult_oi_weights); 
% 
%     % And index of sfGFP
%     species_oi_idx_4lowConc = find(strcmp(fitted_data_4lowConc(1).DataNames,species_name_oi));
%     species_oi_idx_3highConc = find(strcmp(fitted_data_3highConc(1).DataNames,species_name_oi));
%     species_oi_idx_weights = find(strcmp(fitted_data_weights(1).DataNames,species_name_oi));
% 
%     % Plot out weighted fitting, low/high conc fitting, and experimental
%     % data 
%     figure; 
%     for conc_idx = 1:length(conc_vec)
%         subplot(3,3,conc_idx)
% 
%         % Fitted data 
%         plot(fitted_data_weights(conc_idx).Time,fitted_data_weights(conc_idx).Data(:,species_oi_idx_weights),'LineWidth',1.5,'Color','r','LineStyle','--')
%         hold on 
%         plot(fitted_data_weights(conc_idx + length(conc_vec)).Time,fitted_data_weights(conc_idx + length(conc_vec)).Data(:,species_oi_idx_weights),'LineWidth',1.5,'Color','r')
%         if conc_idx <= 4 
%             plot(fitted_data_4lowConc(conc_idx).Time,fitted_data_4lowConc(conc_idx).Data(:,species_oi_idx_4lowConc),'LineWidth',1.5,'Color','g','LineStyle','--')
%             plot(fitted_data_4lowConc(conc_idx + 4).Time,fitted_data_4lowConc(conc_idx + 4).Data(:,species_oi_idx_4lowConc),'LineWidth',1.5,'Color','g')
%         else
%             plot(fitted_data_3highConc(conc_idx - 4).Time,fitted_data_3highConc(conc_idx - 4).Data(:,species_oi_idx_3highConc),'LineWidth',1.5,'Color','b','LineStyle','--')
%             plot(fitted_data_3highConc(conc_idx - 1).Time,fitted_data_3highConc(conc_idx - 1).Data(:,species_oi_idx_3highConc),'LineWidth',1.5,'Color','b')
%         end
% 
%         % Experimental data 
%         group_num_no_empty = 168 + conc_idx; 
%         group_num_empty = 168 + 7 + conc_idx; 
%         no_empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_no_empty); 
%         no_empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_no_empty); 
%         empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_empty); 
%         empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_empty); 
%         plot(no_empty_timeVec,no_empty_GFP_conc,'Color','k','LineWidth',1.5,'LineStyle','--')
%         hold on
%         plot(empty_timeVec,empty_GFP_conc,'Color','k','LineWidth',1.5)
% 
%         title(sprintf('%.1f Reporter Plasmid',conc_vec(conc_idx)))
%     end
% end
% legend('scaled all conc - no empty','scaled all conc - empty','unscaled low/high conc - no empty','unscaled low/high conc - empty','experimental data - no empty','experimental data - empty')
% 
% RFM_ctrl_result_file = load('param_est_run_save/20230709_param_est_run1945_RFM_integration_ctrl.mat'); 

% % Compare parameter distribution 
% num_iter_pos_crosstalk_captured = 4; 
% for idx = 1:length(sort_idx_4lowConc)
%     idx_oi = sort_idx_4lowConc(idx); 
% end


%% Troubleshoot root of duplicate reactions 

% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject_ctrl = setProblemObject_v2(group_description,param_info_path,[]); 
% all_model_rxn_name = {problemObject_ctrl.Model.Reactions.Reaction};
% [all_model_rxn_name_unique,unique_idx_list,~] = unique(all_model_rxn_name); 
% not_unique_idx_list = setdiff(1:length(all_model_rxn_name),unique_idx_list); 
% dup_rxn_name_list = all_model_rxn_name(not_unique_idx_list); 
    % Duplicate reactions are added in the build model phase 

    % % Define components 
    % tube1 = txtl_extract('plasmid_crosstalk');
    % tube2 = txtl_buffer('plasmid_crosstalk');
    % tube3 = txtl_newtube('gene_expression');
    % 
    % % Add DNA template 
    %     % 3WJdB
    % dna_broc_T7_strong = txtl_add_dna(tube3, 'pT7(23)', 'utrbroc(152)', 'no_protein',0, 'plasmid');	
    % dna_broc_T7_weak = txtl_add_dna(tube3, 'pT773(23)', 'utrbroc(152)', 'no_protein',0, 'plasmid');	
    % dna_broc_sigma70_strong = txtl_add_dna(tube3, 'pJ23119(35)', 'utrbroc(800)', 'no_protein',0, 'plasmid');
    %     % sfGFP
    % dna_GFP_T7_strong = txtl_add_dna(tube3,'pT7(23)', 'utrGFP(57)', 'sfGFP(723)',0,'plasmid'); 
    % dna_GFP_T7_weak = txtl_add_dna(tube3,'pT773(23)', 'utrGFP(57)', 'sfGFP(723)',0,'plasmid'); 
    % dna_GFP_sigma70_strong = txtl_add_dna(tube3,'pJ23119(35)', 'utrGFP(57)', 'sfGFP(723)',0,'plasmid'); 
    % dna_GFP_sigma70_weak = txtl_add_dna(tube3,'pJ23105(35)', 'utrGFP(57)', 'sfGFP(723)',0,'plasmid'); 
    %     % empty plasmids 
    % dna_kan = txtl_add_dna(tube3,'pkanR(20)','utrkanR(33)','kanR(783)',0,'plasmid');
    % dna_empty_T7 = txtl_add_dna(tube3,'pT7(23)','utrempty(51)','no_protein',0,'plasmid'); 
    % dna_empty_sigma70 = txtl_add_dna(tube3,'pJ23119(35)','utrempty(51)','no_protein',0,'plasmid'); 
    % 
    % Mobj = txtl_combine([tube1, tube2, tube3]);
    % % txtl_runsim mode initializes reactions and parameters to model object 
    % 
    % [initial_simData] = txtl_runsim(Mobj,14*60*60);
    % 
    % % Add fields to problemObject
    % problemObject.Model = Mobj; 

% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject_ctrl = setProblemObject_v2(group_description,param_info_path,[]); 
% problemObject_updated = setProblemObject_v3(group_description,[],param_info_path); 
% 
% model_rxn_ctrl = problemObject_ctrl.Model.Reactions; 
% model_rxn_updated = problemObject_updated.Model.Reactions; 


%% Test impact of duplicate species and reactions 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject_ctrl = setProblemObject_v2(group_description,param_info_path,[]); 
% problemObject_remove_dup = setProblemObject_v2(group_description,param_info_path,[]); 
% model_remove_dup = problemObject_remove_dup.Model; 
% model_species_name = {model_remove_dup.Species.Name}; 
% model_rxn_name = {model_remove_dup.Reactions.Reaction};
% unique_model_species_name = unique(model_species_name); 
% unique_model_rxn_name = unique(model_rxn_name); 
%     % Check difference in reaction kinetic parameter and species for each
%     % recurring reaction
%         % nan if no duplicate, 0 if not the same, 1 if the same 
% same_rxn_kinetics = nan(length(model_remove_dup.Reactions),1); 
% same_rxn_kineticLaw_Parameters = nan(length(model_remove_dup.Reactions),1); 
% same_rxn_Parameters = nan(length(model_remove_dup.Reactions),1); 
% 
% for rxn_idx = 1:length(model_remove_dup.Reactions)
% 
%     rxn_name = model_remove_dup.Reactions(rxn_idx).Reaction; 
%     rxn_name_idx_list = find(strcmp(rxn_name,model_rxn_name)); 
%     if length(rxn_name_idx_list) > 1 
%         rxn_oi_1 = model_remove_dup.Reactions(rxn_name_idx_list(1)); 
%         rxn_oi_2 = model_remove_dup.Reactions(rxn_name_idx_list(2)); 
%         % Compare fields in the reaction object 
%             % Reaction kinetics 
%         same_rxn_kinetics(rxn_idx) = strcmp(rxn_oi_1.ReactionRate,rxn_oi_2.ReactionRate);
%             % KineticLaw/ Parameters 
%         rxn_oi_1_kineticLaw = rxn_oi_1.KineticLaw; 
%         rxn_oi_2_kineticLaw = rxn_oi_2.KineticLaw; 
%         rxn_oi_1_kineticParameters = rxn_oi_1_kineticLaw.ParameterVariableNames; 
%         rxn_oi_2_kineticParameters = rxn_oi_2_kineticLaw.ParameterVariableNames; 
%         rxn_oi_1_kineticLawParams = {rxn_oi_1_kineticLaw.Parameters.Name}; 
%         rxn_oi_2_kineticLawParams = {rxn_oi_2_kineticLaw.Parameters.Name}; 
%         same_rxn_Parameters(rxn_idx) = all(strcmp(rxn_oi_1_kineticParameters,rxn_oi_2_kineticParameters));
%         same_rxn_kineticLaw_Parameters(rxn_idx) = all(strcmp(rxn_oi_1_kineticLawParams,rxn_oi_2_kineticLawParams)); 
%     end
% 
% end

% All duplicate reactions are same in terms of reaction kinetics and
% kinetic parameters 


%     % Delete redundant reactions if all kinetic laws and parameters are the
%     % same. Let's delete the second appearing one 
% rxn_idx = 1;
% while rxn_idx < length(model_remove_dup.Reactions)
%     rxn_oi = model_remove_dup.Reactions(rxn_idx); 
%     rxn_name = rxn_oi.Reaction; 
%     all_model_rxn_name = {model_remove_dup.Reactions.Reaction}; 
%     rxn_name_idx_list = find(strcmp(rxn_name,all_model_rxn_name)); 
%     if length(rxn_name_idx_list) > 1 
%         rxn_to_be_deleted = model_remove_dup.Reactions(rxn_name_idx_list(2)); 
%         delete(rxn_to_be_deleted); 
%     end
% 
%     rxn_idx = rxn_idx + 1; 
% end
% 
%     % Create a simFunction and dosing information
% observable_species_name = 'protein sfGFP*'; 
% observable_species_idx = find(strcmp(observable_species_name,model_species_name)); 
% additional_track_species = [model_species_name(1:observable_species_idx - 1),model_species_name(observable_species_idx + 1:end)]; 
% simFunction_ctrl = create_simFun_from_problemObject(problemObject_ctrl,additional_track_species);
% simFunction_remove_dup = create_simFun_from_problemObject(problemObject_remove_dup,additional_track_species); 
% dosing_information = create_dosing_info_from_problemObject(problemObject_ctrl); 
% 
% %     % Run a parameter estimation control 
% % [fitResults_ctrl,simdata_ctrl] = fit(problemObject_ctrl);
% % [fitResults_remove_dup,simdata_remove_dup] = fit(problemObject_remove_dup);
% 
%     % Run simFunction using initial parameter values 
% estimated_params = [problemObject_ctrl.Estimated.InitialValue]; 
% [simulated_time_ctrl,simulated_data_ctrl] = simFunction_ctrl(estimated_params,tEnd,dosing_information,tStart:tEnd);
% [simulated_time_remove_dup,simulated_data_remove_dup] = simFunction_remove_dup(estimated_params,tEnd,dosing_information,tStart:tEnd);
% 
%     % Compare simulated data 
%         % Intrapolate them to nT = 100 
% same_time_course_flag = nan(length(simulated_time_ctrl),1); 
% allowed_deviation = 1e-04; 
% target_time_interval = (tEnd - tStart) / 1000; 
% target_timeVec = tStart:target_time_interval:tEnd; 
% for simData_idx = 1:length(simulated_time_ctrl)
% 
%     % Select time vec and time course 
%     simData_time_ctrl_single = simulated_time_ctrl{simData_idx}; 
%     simData_time_remove_dup_single = simulated_time_remove_dup{simData_idx}; 
%     simData_data_ctrl_single = simulated_data_ctrl{simData_idx};
%     simData_data_remove_dup_single = simulated_data_remove_dup{simData_idx}; 
% 
%     % Remove duplicate values for intrapolation
%     [simData_time_ctrl_single_unique,ctrl_unique_idx,~] = unique(simData_time_ctrl_single); 
%     [simData_time_remove_dup_single_unique,remove_dup_unique_idx,~] = unique(simData_time_remove_dup_single); 
%     simData_data_ctrl_single_unique = simData_data_ctrl_single(ctrl_unique_idx,:); 
%     simData_data_remove_dup_single_unique = simData_data_remove_dup_single(remove_dup_unique_idx,:); 
% 
%     % Intrapolate to target timeVec
%     simData_data_ctrl_single_interp = interp1(simData_time_ctrl_single_unique,simData_data_ctrl_single_unique,target_timeVec,"linear",'extrap');
%     simData_data_remove_dup_single_interp = interp1(simData_time_remove_dup_single_unique,simData_data_remove_dup_single_unique,target_timeVec,"linear",'extrap');
% 
%     % Compare 
%     same_time_course_flag(simData_idx) = all(all(abs(simData_data_remove_dup_single_interp - simData_data_ctrl_single_interp) < allowed_deviation)); 
% 
% 
% end
% 
% 
% figure; 
% plot(target_timeVec,simData_data_ctrl_single_interp(:,1),'k')
% hold on
% plot(target_timeVec,simData_data_remove_dup_single_interp(:,1),'r')
%% Test RFM integration results 
% RFM_integration_ctrl_result_file = load('param_est_run_save/20230704_param_est_run0750_RFM_integration_ctrl.mat'); 
% RFM_integration_test_result_file = load('param_est_run_save/20230703_param_est_run1546_RFM_integration.mat'); 

%     % Compile result file 
% wrapper_analyze_fit_result_v2(RFM_integration_ctrl_result_file.all_fitResults,RFM_integration_ctrl_result_file.problemObject,'param_est_run_save/20230704_param_est_run0750_RFM_integration_ctrl.mat','SSE',false,false) 
% wrapper_analyze_fit_result_v2(RFM_integration_test_result_file.all_fitResults,RFM_integration_test_result_file.problemObject,'param_est_run_save/20230703_param_est_run1546_RFM_integration.mat','SSE',false,false) 

% Check fitted SSE 
% figure;
% histogram(log10(RFM_integration_ctrl_result_file.metric_summary));
% figure; 
% histogram(log10(RFM_integration_test_result_file.metric_summary)); 
% 
% % Plot out the experimental data 
% problemObject_ctrl = RFM_integration_ctrl_result_file.problemObject;
% % problemObject_test = RFM_integration_test_result_file.problemObject; 
% experimental_groupedData = problemObject_ctrl.Data; 
% figure; 
% for conc_idx = 1:length(conc_vec)
%     subplot(3,3,conc_idx)
%     group_num_no_empty = 168 + conc_idx; 
%     group_num_empty = 168 + 7 + conc_idx; 
%     no_empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_no_empty); 
%     no_empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_no_empty); 
%     empty_timeVec = experimental_groupedData.Time(experimental_groupedData.Group == group_num_empty); 
%     empty_GFP_conc = experimental_groupedData.GFP_concentration(experimental_groupedData.Group == group_num_empty); 
%     plot(no_empty_timeVec,no_empty_GFP_conc,'Color','g','LineWidth',1.5)
%     hold on
%     plot(empty_timeVec,empty_GFP_conc,'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     title(num2str(conc_vec(conc_idx))); 
% 
% end
% 
% 
% % For each result file, select top 5 best fitted model and check crosstalk
% % ratio 
% [~,sorted_SSE_idx_ctrl] = sort(RFM_integration_ctrl_result_file.metric_summary,'ascend'); 
% [~,sorted_SSE_idx_test] = sort(RFM_integration_test_result_file.metric_summary,'ascend'); 
% 
% species_oi = 'protein sfGFP*'; 
% 
% for idx = 1:5
%     % Select fitResult of interest and get simData
%     fitResult_oi_ctrl = RFM_integration_ctrl_result_file.all_fitResults{sorted_SSE_idx_ctrl(idx)}; 
%     fitResult_oi_test = RFM_integration_test_result_file.all_fitResults{sorted_SSE_idx_test(idx)}; 
%     simData_ctrl = fitted(fitResult_oi_ctrl);
%     simData_test = fitted(fitResult_oi_test); 
% 
%     % These should be organized as no empty - conc_vec 1-7, empty -
%     % conc_vec 1-7 ? Let's try this first and plot out the time-course for
%     % GFP 
%     figure;
%     for conc_idx = 1:length(conc_vec)
% 
%         subplot(3,3,conc_idx)
%         simData_no_empty_ctrl = simData_ctrl(conc_idx); 
%         simData_no_empty_test = simData_test(conc_idx);
%         simData_empty_ctrl = simData_ctrl(length(conc_vec) + conc_idx); 
%         simData_empty_test = simData_test(length(conc_vec) + conc_idx);
% 
%         species_oi_idx_ctrl = find(strcmp(simData_no_empty_ctrl.DataNames,species_oi));
%         species_oi_idx_test = find(strcmp(simData_no_empty_test.DataNames,species_oi));
% 
% 
%         plot(simData_no_empty_ctrl.Time,simData_no_empty_ctrl.Data(:,species_oi_idx_ctrl),'--','LineWidth',1.5,'Color','g')
%         hold on 
%         plot(simData_no_empty_test.Time,simData_no_empty_test.Data(:,species_oi_idx_test),'-','LineWidth',1.5,'Color','g')
%         plot(simData_empty_ctrl.Time,simData_empty_ctrl.Data(:,species_oi_idx_ctrl),'--','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         plot(simData_empty_test.Time,simData_empty_test.Data(:,species_oi_idx_test),'-','LineWidth',1.5,'Color',[0.5,0.5,0.5])
% 
%     end
%     legend('Ctrl - no empty','Test - no empty','Ctrl - empty','Test - empty')
% end

% % Check estimated parameter distribution 
% num_params_ctrl = height(RFM_integration_ctrl_result_file.all_fitResults{1}.ParameterEstimates); 
% num_params_test = height(RFM_integration_test_result_file.all_fitResults{1}.ParameterEstimates); 
% 
%     % Compile estimated parameters 
% all_estimated_params_ctrl = nan(num_iter,num_params_ctrl); 
% all_estimated_params_test = nan(num_iter,num_params_test); 
% 
% for iter = 1:num_iter
% 
%     estimated_params_ctrl = [RFM_integration_ctrl_result_file.all_fitResults{iter}.ParameterEstimates.Estimate];
%     estimated_params_test = [RFM_integration_test_result_file.all_fitResults{iter}.ParameterEstimates.Estimate];
% 
%     all_estimated_params_ctrl(iter,:) = estimated_params_ctrl;
%     all_estimated_params_test(iter,:) = estimated_params_test; 
% 
%     if isequal(iter,1)
%         estimated_params_name_ctrl = RFM_integration_ctrl_result_file.all_fitResults{iter}.ParameterEstimates.Name;
%         estimated_params_name_test = RFM_integration_test_result_file.all_fitResults{iter}.ParameterEstimates.Name;
%     end
% 
% end

% param_name_oi_list = {'TXTL_UTR_kanR_R','TXTL_UTR_GFP_R'};
% param_name_oi_idx_list = [10,17];
% 
% all_UTR_kanR_R_ctrl = all_estimated_params_ctrl(:,17); 
% all_UTR_GFP_R_ctrl = all_estimated_params_ctrl(:,10); 
% all_UTR_kanR_R_test = all_estimated_params_test(:,17); 
% all_UTR_GFP_R_test = all_estimated_params_test(:,10); 

%     % Check parameter distribution
% figure;
% for param_idx = 1:length(estimated_params_name_test)
%     estimated_param_name = estimated_params_name_test{param_idx};
%     estimated_param_test = all_estimated_params_test(:,param_idx); 
%     subplot(ceil(sqrt(length(estimated_params_name_test))),ceil(sqrt(length(estimated_params_name_test))),param_idx)
%     histogram(log10(estimated_param_test)); 
%     hold on 
%     if ~isequal(param_idx,length(estimated_params_name_test))
%         estimated_param_ctrl = all_estimated_params_ctrl(:,param_idx); 
%         histogram(log10(estimated_param_ctrl)); 
%     end
%     title(strrep(estimated_param_name,'_',' '));
% 
% end

% % Check ribosome usage in both cases - if there's no increased competition for
% % ribosomes the added mechanism won't work 
% [~,sorted_SSE_idx_ctrl] = sort(RFM_integration_ctrl_result_file.metric_summary,'ascend'); 
% [~,sorted_SSE_idx_test] = sort(RFM_integration_test_result_file.metric_summary,'ascend'); 
% 
% % Define species we're interested in & get the indices
% species_oi_list = {'Ribo','RNase','RNAP','AGTP','CUTP','RNase','toxin'};
%     % Get DataNames 
% sample_fitResult_ctrl = RFM_integration_ctrl_result_file.all_fitResults{1}; 
% sample_fitResult_test = RFM_integration_test_result_file.all_fitResults{1}; 
% sample_fitted_data_ctrl = fitted(sample_fitResult_ctrl);
% sample_fitted_data_test = fitted(sample_fitResult_test); 
% sample_dataNames_ctrl = sample_fitted_data_ctrl(1).DataNames; 
% sample_dataNames_test = sample_fitted_data_test(1).DataNames; 
%     % Find species of interest in DataNames and gather indices 
% species_oi_idx_list_ctrl = nan(length(species_oi_list),1); 
% species_oi_idx_list_test = nan(length(species_oi_list),1); 
% for species_oi_idx = 1:length(species_oi_list)
%     species_oi = species_oi_list{species_oi_idx}; 
%     species_oi_idx_list_ctrl(species_oi_idx) = find(strcmp(sample_dataNames_ctrl,species_oi)); 
%     species_oi_idx_list_test(species_oi_idx) = find(strcmp(sample_dataNames_test,species_oi)); 
% end
% 
% for idx = 1:5
%     % Get fitted data 
%     fitResult_oi_ctrl = RFM_integration_ctrl_result_file.all_fitResults{idx};
%     fitResult_oi_test = RFM_integration_test_result_file.all_fitResults{idx}; 
%     fitted_data_ctrl = fitted(fitResult_oi_ctrl);
%     fitted_data_test = fitted(fitResult_oi_test); 
%     fitted_data_ctrl_empty = fitted_data_ctrl(end);
%     fitted_data_test_empty = fitted_data_test(end); 
%     fitted_data
% 
%     figure; 
%     gen
%     for species_idx = 1:length(species_oi_list)
%         subplot(ceil(sqrt(length(species_oi_list))),ceil(sqrt(length(species_oi_list))),species_idx)
%         plot(fitted_data_ctrl_empty.Time,fitted_data_ctrl_empty.Data(:,species_oi_idx_list_ctrl(species_idx)),'LineWidth',1.5,'Color','k');
%         hold on 
%         plot(fitted_data_test_empty.Time,fitted_data_test_empty.Data(:,species_oi_idx_list_test(species_idx)),'LineWidth',1.5,'Color','r');
%         title(species_oi_list(species_idx))
%     end
% 
% end

%% Rename previous sampling results 
% for iter = 21:1000
%     result_fileName = sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d.mat',iter); 
%     if exist(result_fileName,'file')
%         new_result_fileName = sprintf('param_est_run_save/20240701_param_and_init_cond_sampling_constrained_run%d_old.mat',iter);
%         movefile(result_fi    leName,new_result_fileName); 
%     end
% end


%% In previous parameter estimation attempts, how do the TXTL_UTR's compare to each other? 
% large_scale_result_file = load('param_est_run_save/20231106_large_scale_run.mat'); 
% sample_result_file = load('20231102_param_est_run1534_100.mat');
% 
% param_name_list = {sample_result_file.problemObject.Estimated.Name}; 
% param_name_oi_list = {'TXTL_UTR_kanR_R','TXTL_UTR_GFP_R'};
% param_name_oi_idx_list = nan(length(param_name_oi_list),1); 
%     % Get the parameter index in list 
% for param_idx = 1:length(param_name_oi_list)
%     param_name_oi = param_name_oi_list{param_idx}; 
%     param_name_oi_idx_list(param_idx,1) = find(strcmp(param_name_oi,param_name_list)); 
% end
% 
%     % Get the parameters out from the matrix 
% param_val_oi = nan(height(large_scale_result_file.updated_optimization_run_table),length(param_name_oi_idx_list)); 
% all_estimated_params = large_scale_result_file.updated_optimization_run_table.ParameterEstimates; 
% for iter_idx = 1:height(large_scale_result_file.updated_optimization_run_table)
%     estimated_params = all_estimated_params{iter_idx,1};
%     param_val_oi(iter_idx,:) = estimated_params(param_name_oi_idx_list); 
% 
% end
% 
% % Let's visualize these by plotting out the ratios? Just keep in mind that
% % it will be skewed to the larger than 1 side? 
%     % Expect TXTL_UTR_kanR_R > TXTL_UTR_GFP_R (stronger binding) 
% % ratio_kanR_R_over_GFP_R = param_val_oi(:,1) ./ param_val_oi(:,2); % Expect to be larger than 1 
% % histogram(ratio_kanR_R_over_GFP_R)
%     % Let's do bar plot of # samples vs. percentage difference? 
% percentage_diff = (param_val_oi(:,2) - param_val_oi(:,1)) ./ param_val_oi(:,2); 
% edges = [-100,-0.5,-0.25,-0.1,0.1,0.25,0.5,100];
% binned_percentage_diff = discretize(percentage_diff,edges); 
% num_in_category = nan(length(edges)-1,1); 
% for edge_idx = 1:length(edges)-1
%     num_in_category(edge_idx) = length(find(binned_percentage_diff == edge_idx)); 
% end
% Xlabels = {'>50% Smaller','>25% Smaller','>10% Smaller','Within 10% deviation','>10% Larger','>25% Larger','>50% Larger'}; 
% X = categorical(Xlabels);
% X = reordercats(X,Xlabels);
% bar(X,num_in_category)


%% Check sampling saved files to check what can be dropped out to optimiza storage 
% test_result_file = load('test_save_files/20240703_temp_check_storage.mat'); 
% sampled_params_selected = test_result_file.sampled_params_selected;
% all_dosing_information = test_result_file.all_dosing_information; 
% all_simulated_time = test_result_file.all_simulated_time;
% all_simulated_data = test_result_file.all_simulated_data; 
% penalty_values = test_result_file.penalty_values; 
% penalty_labels = test_result_file.penalty_labels; 
% save('test_save_files/20240703_temp_check_save2.mat','sampled_params_selected','all_dosing_information','all_simulated_time','all_simulated_data',...
%     'penalty_values',"penalty_labels");

%% Sequential toxin mechanism fitting 

% group_description = {'PE_T7_strong_no_empty','PE_T7_weak_no_empty','PE_sigma70_strong_no_empty','PE_sigma70_weak_no_empty'};
% no_toxin_group_number = [1:6,29:34,57:61,85:91];
% toxin_group_number = [7,35,62,63]; 
% param_info_path = 'parameters_v2.xlsx';
% 
% % Fitting all no toxin impact datasets 
% problemObject = setProblemObject_v2_no_toxin(group_description,no_toxin_group_number,param_info_path);
% save_path_id = 'sequential_toxin'; 
% wrapper_fit_to_data_w_probObject(problemObject,save_path_id); 
% 
% % Fit individual dose-response curve 
% T7_strong_group_description = {'PE_T7_strong_no_empty'}; 
% T7_weak_group_description = {'PE_T7_weak_no_empty'};
% sigma70_strong_group_description = {'PE_sigma70_strong_no_empty'};
% sigma70_weak_group_description = {'PE_sigma70_weak_no_empty'};
% 
% % problemObject_T7_strong = setProblemObject_v2(T7_strong_group_description,param_info_path,[]);
% % problemObject_T7_weak = setProblemObject_v2(T7_weak_group_description,param_info_path,[]);
% % problemObject_sigma70_strong = setProblemObject_v2(sigma70_strong_group_description,param_info_path,[]);
% problemObject_sigma70_weak = setProblemObject_v2(sigma70_weak_group_description,param_info_path,[]);
% problemObject_all_baseline = setProblemObject_v2(group_description,param_info_path,[]);
% 
% wrapper_fit_to_data_w_probObject(problemObject_all_baseline,'all_baseline'); 

% sequential_toxin_result_file_name = 'param_est_run_save/20230616_param_est_run1130_sequential_toxin.mat';
% T7_strong_toxin_result_file_name = 'param_est_run_save/20230614_param_est_run0906_baseline_T7_strong.mat'; 
% T7_weak_toxin_result_file_name = 'param_est_run_save/20230614_param_est_run1027_baseline_T7_weak.mat'; 
% sigma70_weak_result_file_name = 'param_est_run_save/20230616_param_est_run1258_baseline_sigma70_weak.mat';
% 
% sequential_toxin_result_file = load(sequential_toxin_result_file_name); 
% T7_strong_toxin_result_file = load(T7_strong_toxin_result_file_name);
% T7_weak_toxin_result_file = load(T7_weak_toxin_result_file_name);
% sigma70_weak_toxin_result_file = load(sigma70_weak_result_file_name); 

% % Check fitting results and compare parameter distributions
% all_fitResults_T7_strong = T7_strong_toxin_result_file.all_fitResults;
% all_fitResults_T7_weak = T7_weak_toxin_result_file.all_fitResults;
% all_fitResults_sigma70_weak = sigma70_weak_result_file.all_fitResults;
% all_fitResults_sequential_toxin = sequential_toxin_result_file.all_fitResults; 
% 
% [all_estimated_params_T7_strong,all_SSE_T7_strong] = compile_fitResult(all_fitResults_T7_strong); 
% [all_estimated_params_T7_weak,all_SSE_T7_weak] = compile_fitResult(all_fitResults_T7_weak); 
% [all_estimated_params_sigma70_weak,all_SSE_sigma70_weak] = compile_fitResult(all_fitResults_sigma70_weak); 
% [all_estimated_params_sequential_toxin,all_SSE_sequential_toxin] = compile_fitResult(all_fitResults_sequential_toxin); 

% figure; 
% histogram(log10(all_SSE_T7_strong))
% title('T7 strong SSE distribution')
% plot_param_distribution({T7_strong_toxin_result_file.problemObject.Estimated.Name},all_estimated_params_T7_strong,'T7 strong parameter distribution')
% 
% figure; 
% histogram(log10(all_SSE_T7_weak))
% title('T7 weak SSE distribution')
% plot_param_distribution({T7_weak_toxin_result_file.problemObject.Estimated.Name},all_estimated_params_T7_weak,'T7 weak parameter distribution')
% 
% figure; 
% histogram(log10(all_SSE_sigma70_weak))
% title('sigma70 weak SSE distribution')
% plot_param_distribution({sigma70_weak_result_file.problemObject.Estimated.Name},all_estimated_params_sigma70_weak,'sigma70 weak parameter distribution')
% 
% figure; 
% histogram(log10(all_SSE_sequential_toxin))
% title('sequential toxin SSE distribution')
% plot_param_distribution(all_SSE_sequential_toxin,{sequential_toxin_result_file.problemObject.Estimated.Name},all_estimated_params_sequential_toxin,'sequential toxin parameter distribution')

% % Compare fitting 
% [~,sort_idx] = sort(all_SSE_T7_strong,'ascend'); 
% for idx = 1:5
%     fitResult_oi = all_fitResults_T7_strong{sort_idx(idx)};
%     simData = fitted(fitResult_oi); 
%     subplot(3,3,idx)
%     for conc_idx = 1:length(simData)
%         simData_single = simData(conc_idx);
%         plot(simData_single.Time,simData_single.Data(:,data_idx))
%         hold on 
% 
%     end
% end


% function plot_param_distribution(all_SSE,param_names,param_values,plot_title)
%     figure;
%     selected_iter_idx_list = all_SSE < 1e+12; 
%     unselected_iter_idx_list = all_SSE > 1e+12;
%     for param_idx = 1:length(param_names)
%         param_name = param_names{param_idx};
%         selected_param_vals = param_values(selected_iter_idx_list,param_idx); 
%         histogram(selected_param_vals)
%         hold on 
%         unselected_param_vals = param_values(unselected_iter_idx_list,param_idx); 
%         histogram(unselected_param_vals)
%         subplot(ceil(sqrt(length(param_names))),ceil(sqrt(length(param_names))),param_idx)
%         title(strrep(param_name,'_',' '))
%     end
%     sgtitle(plot_title)
% end
% 
% function [all_estimated_params,all_SSE] = compile_fitResult(all_fitResults)
% 
%     fitResult_test = all_fitResults{1}; 
%     num_params = height(fitResult_test.ParameterEstimates); 
%     num_iter = length(all_fitResults);
% 
%     all_estimated_params = nan(num_iter,num_params); 
%     all_SSE = nan(num_iter,1); 
%     for iter = 1:length(all_fitResults)
%         fitResult_oi = all_fitResults{iter}; 
%         all_SSE(iter) = fitResult_oi.SSE;
%         all_estimated_params(iter,:) = [fitResult_oi.ParameterEstimates.Estimate];
%     end
% 
% end
%% Check partial fitting results + parameter distribution
% all_modes = {'TX','PE'};
% all_combination_list = {'no_empty','empty','empty_T7','empty_sigma70'}; 
% for mode_idx = 2:length(all_modes)
%     mode = all_modes{mode_idx};
%     for promotor_idx = 1:length(promotor_name_list)
%         promotor_name = promotor_name_list{promotor_idx};
%         for combo_idx = 1:length(all_combination_list)
%             combo_name = all_combination_list{combo_idx};
%             result_file_name = sprintf('param_est_run_save/20230813_param_est_run0807_%s_%s_%s.mat',...
%                 mode,promotor_name,combo_name); 
%             result_file = load(result_file_name);
%             % Figure 1 - SSE distribution 
%             figure; 
%             histogram(log10(result_file.metric_summary))
%             xlabel('log_{10}SSE')
%             ylabel('# Occurence')
%             title(strrep(sprintf('%s Mode %s %s',mode,promotor_name,combo_name),'_',' '))
% 
%             % Figure 2 - parameter distribution 
%             figure;
% 
%             param_name_list = result_file.estimated_param_summary.Name;
%             all_estimated_params_mat = nan(num_iter,length(param_name_list));
%             for iter = 1:num_iter
%                 eval(sprintf('fitted_params = [result_file.estimated_param_summary.Estimate_%d];',iter)); 
%                 all_estimated_params_mat(iter,:) = fitted_params';
%             end
% 
%             for param_idx = 1:length(param_name_list)
%                 param_name = param_name_list{param_idx}; 
%                 subplot(ceil(sqrt(length(param_name_list))),ceil(sqrt(length(param_name_list))),param_idx)
%                 histogram(log10(all_estimated_params_mat(:,param_idx))); 
%                 title(strrep(param_name,'_',' '))
%                 xlabel('log_{10}parameter value')
%                 ylabel('# Occurence')
%             end
%             sgtitle(strrep(sprintf('%s Mode %s %s',mode,promotor_name,combo_name),'_',' '))
% 
% 
%         end
%     end
% end

%% Check goal attainment results /objString = 20,22 results 

% % Get all file names 
% directory_path = 'param_est_run_save'; 
% directory_files = dir(directory_path); 
% file_names = {directory_files.name};
% 
% num_iter = 48; 
% 
%     % Goal attainment 
% % all_fitResult = cell(num_iter,1); 
% % all_attain_factor = nan(num_iter,1); 
% % all_scaled_init_penalty_terms = nan(num_iter,33); 
% % all_scaled_final_penalty_terms = nan(num_iter,33); 
% % all_estimated_params = nan(num_iter,31); 
% % all_goal_attainment = nan(num_iter,33); 
% 
%     % single objective objIdx = 20 
% all_fitResult = cell(num_iter,1); 
% all_fval = nan(num_iter,1); 
% all_estimated_params = nan(num_iter,31); 
% all_penalty_init = nan(num_iter,7);
% all_penalty_final = nan(num_iter,7); 
% 
% for iter = 1:num_iter
%     % resultFileName_idx = contains(file_names,sprintf('goal_attainment_run1_iter_%d',iter)); 
%     resultFileName_idx = contains(file_names,sprintf('objString_20_iter_%d',iter)); 
%     if any(resultFileName_idx)
%         resultFileName = file_names{resultFileName_idx};
%         resultFile = load(resultFileName); 
% 
%         if ~exist('simFunction','var')
%             problemObject = resultFile.problemObject;
%             dosing_information = create_dosing_info_from_problemObject(problemObject); 
%             simFunction = create_simFun_from_problemObject(problemObject);
%         end
%         % all_fitResult{iter,1} = resultFile.fitResult;
%         % all_attain_factor(iter,1) = resultFile.attainfactor;
%         % all_scaled_init_penalty_terms(iter,:) = resultFile.scaled_penalty_vec_init';
%         % all_scaled_final_penalty_terms(iter,:) = resultFile.scaled_penalty_vec_final';
%         % all_estimated_params(iter,:) = resultFile.estimated_params; 
%         % goal_attainment = resultFile.fval - resultFile.goal; % Satisfied = < 0
%         % all_goal_attainment(iter,:) = goal_attainment'; 
% 
%         all_fitResult{iter,1} = resultFile.fitResult; 
%         all_fval(iter,1) = resultFile.fitResult.fval; 
%         all_estimated_params(iter,:) = resultFile.fitResult.estimated_params; 
%         all_penalty_init(iter,:) = resultFile.penalty_term_init; 
%         all_penalty_final(iter,:) = resultFile.penalty_term_final; 
% 
%         % estimated_params = resultFile.estimated_params; 
%         % [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
%         % plot_simulated(simulated_time,simulated_data,'baseline')
%         % plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
%     end
% 
% end
% 
% % Plot out the objective value 
% figure;
% histogram(all_fval);
% xlabel('Large positive crosstalk penalty')
% ylabel('# Occurence')
% title('Large Positive Crosstalk Penalty Distribution')


% % % Plot out the difference between max weak promotor positive crosstalk
% % % ratio and max strong promotor positive crosstalk ratio 
% % figure;
% % nbins = 20; 
% % h = histogram(all_fval,nbins);
% % xlabel('max(T7, sigma70 strong positivie crosstalk ratio) - max(sigma70 weak positive crosstalk ratio)'); 
% % ylabel('# Occurence')
% % title('Optimizing for large positive crosstalk ratio')
% 
% % Check out crosstalk ratio plots for the most negative ones 
% 
% [sorted_fval,sort_idx_list] = sort(all_fval,'ascend'); 
% for k = 1:10
% 
%     fitResult_oi = all_fitResult{sort_idx_list(k)};
%     estimated_params = all_estimated_params(sort_idx_list(k),:); 
%     [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% 
% end

% % Plot change in individual penalty term throughout the 48 iterations 
% figure; 
% penalty_term_labels = resultFile.penalty_term_labels; 
% penalty_term_length = resultFile.penalty_term_length;
% 
% all_penalty_for_plot = zeros(3 * num_iter,sum(penalty_term_length));
% X_labels_for_plot = cell(3 * num_iter,1); 
% for iter = 1:num_iter
%     penalty_term_init = all_scaled_init_penalty_terms(iter,:);
%     penalty_term_final = all_scaled_final_penalty_terms(iter,:); 
%     all_penalty_for_plot((iter - 1) * 3 + 1,:) = penalty_term_init; 
%     all_penalty_for_plot((iter - 1) * 3 + 2,:) = penalty_term_final';
%     X_labels_for_plot{(iter - 1) * 3 + 1} = sprintf('#%d init',iter); 
%     X_labels_for_plot{(iter - 1) * 3 + 2} = sprintf('#%d final',iter);
%     X_labels_for_plot{(iter - 1) * 3 + 3} = num2str(iter);  
% end
% X_labels = categorical(X_labels_for_plot); 
% X_labels = reordercats(X_labels,X_labels_for_plot); 

%     % Modify the legends 
% penalty_figure_legends = cell(sum(penalty_term_length),1); 
% penalty_idx = 1;
% actual_idx = 1; 
% while penalty_idx <= length(penalty_term_length)
%     penalty_term_label_single = penalty_term_labels{penalty_idx}; 
%     penalty_term_length_single = penalty_term_length(penalty_idx); 
%     if isequal(penalty_term_length_single,1)
%         penalty_figure_legends{actual_idx} = penalty_term_label_single;
%         actual_idx = actual_idx + 1; 
%     else
%         for j = 1:penalty_term_length_single
%             penalty_figure_legends{actual_idx} = sprintf('%s #%d',penalty_term_label_single,j); 
%             actual_idx = actual_idx + 1; 
%         end
%     end
%     penalty_idx = penalty_idx + 1; 
% end
% 
% bar(X_labels,all_penalty_for_plot,'stacked')
% legend(strrep(penalty_figure_legends,'_',' '))
% sgtitle(sprintf('objective string %d',obj_string_idx))

%% Penalty distribution for large-scale fitting
% large_scale_result_file_name = 'test_save_files/20240531_large_scale_fitting_penalty_term.mat';
% large_scale_result_file = load(large_scale_result_file_name);
% 
%     % Check parameter distribution in separate peaks 
%     % Get basic information on the result struct 
% sample_estimated_params = cell2mat(large_scale_result_file.large_scale_fitting_result_struct{1,'ParameterEstimates'});
% num_params = length(sample_estimated_params); 
% num_runs = height(large_scale_result_file.large_scale_fitting_result_struct); 
% param_names = large_scale_result_file.estimated_param_names; 
% 
%     % Convert estimated params into the correct size 
% all_estimated_params_cell = [large_scale_result_file.large_scale_fitting_result_struct.ParameterEstimates];
% all_estimated_params_vec = cell2mat(all_estimated_params_cell);
% all_estimated_params_mat = reshape(all_estimated_params_vec,[num_params,num_runs]);
% all_estimated_params_mat = all_estimated_params_mat'; 
% 
%     % Baseline penalty
% baseline_penalty = large_scale_result_file.large_scale_fitting_result_struct.all_log_baseline_data_penalty;
% baseline_penalty_satisfied_idx_list = find(baseline_penalty > 0); 
% baseline_penalty_unsatisfied_idx_list = find(baseline_penalty < 0); 
%         % Get the estimated params
% baseline_penalty_satisfied_estimated_params = all_estimated_params_mat(baseline_penalty_satisfied_idx_list,:); 
% baseline_penalty_unsatisfied_estimated_params = all_estimated_params_mat(baseline_penalty_unsatisfied_idx_list,:); 
%         % Plot out the histogram
% figure; 
% for param_idx = 1:num_params
%     subplot(ceil(sqrt(num_params)),ceil(sqrt(num_params)),param_idx)
% 
%     estimated_param_name = param_names{param_idx}; 
%     estimated_param_name_for_plot = strrep(estimated_param_name,'_',' ');
%     estimated_param_name_for_plot = strrep(estimated_param_name_for_plot,'TXTL',''); 
% 
%         % log transform each parameter? 
%     log_satisfied_estimated_params_vec = log10(baseline_penalty_satisfied_estimated_params(:,param_idx)); 
%     h_satisfied = histogram(log_satisfied_estimated_params_vec);
%     hold on 
% 
%     log_unsatisfied_estimated_params_vec = log10(baseline_penalty_unsatisfied_estimated_params(:,param_idx)); 
%     h_unsatisfied = histogram(log_unsatisfied_estimated_params_vec);
% 
%     title(estimated_param_name_for_plot)
% 
% end
% legend('Satisfied','Unsatisfied')
% sgtitle('Parameter distribution for baseline penalty')
% 
%     % Crosstalk transition penalty 
% crosstalk_transition_penalty =  large_scale_result_file.large_scale_fitting_result_struct.all_crosstalk_ratio_transition_penalty;
%         % Note: Does not mean the transition is present in every case
% crosstalk_transition_satisfied_idx_list = find(crosstalk_transition_penalty < 0);
% crosstalk_transition_unsatisfied_idx_list = find(crosstalk_transition_penalty > 0);
% 
%         % Get the estimated params
% crosstalk_transition_penalty_satisfied_estimated_params = all_estimated_params_mat(crosstalk_transition_satisfied_idx_list,:); 
% crosstalk_transition_penalty_unsatisfied_estimated_params = all_estimated_params_mat(crosstalk_transition_unsatisfied_idx_list,:); 
%         % Plot out the histogram
% figure; 
% for param_idx = 1:num_params
%     subplot(ceil(sqrt(num_params)),ceil(sqrt(num_params)),param_idx)
% 
%     estimated_param_name = param_names{param_idx}; 
%     estimated_param_name_for_plot = strrep(estimated_param_name,'_',' ');
%     estimated_param_name_for_plot = strrep(estimated_param_name_for_plot,'TXTL',''); 
% 
%         % log transform each parameter? 
%     log_satisfied_estimated_params_vec = log10(crosstalk_transition_penalty_satisfied_estimated_params(:,param_idx)); 
%     h_satisfied = histogram(log_satisfied_estimated_params_vec);
%     hold on 
% 
%     log_unsatisfied_estimated_params_vec = log10(crosstalk_transition_penalty_unsatisfied_estimated_params(:,param_idx)); 
%     h_unsatisfied = histogram(log_unsatisfied_estimated_params_vec);
% 
%     title(estimated_param_name_for_plot)
% 
% end
% legend('Satisfied','Unsatisfied')
% sgtitle('Parameter distribution for crosstalk transition penalty')

%     % log(SSE_baseline) 
% figure;
% histogram(log10(large_scale_result_file.large_scale_fitting_result_struct.all_baseline_data_dev),'BinLimits',[10.5,15])
% title('Baseline Deviation')
% xlabel('log_{10}SSE')
% ylabel('# Occurence')
% set(gca,'FontSize',12)
% 
%     % crosstalk ratio 
% figure;
% histogram(large_scale_result_file.large_scale_fitting_result_struct.all_crosstalk_ratio_dev)
% title('Crosstalk Ratio Deviation')
% xlabel('SSE')
% ylabel('# Occurence')
% set(gca,'FontSize',12)
% 
%     % log baseline penalty 
% figure;
% histogram(large_scale_result_file.large_scale_fitting_result_struct.all_log_baseline_data_penalty)
% title('Baseline Data Penalty')
% xlabel('Log baseline data penalty')
% ylabel('# Occurence')
% set(gca,'FontSize',12)
% 
%     % Crosstalk ratio transition penalty 
% figure;
% histogram(large_scale_result_file.large_scale_fitting_result_struct.all_crosstalk_ratio_transition_penalty)
% title('Crosstalk Transition Penalty')
% xlabel('Crosstalk Transition Penalty')
% ylabel('# Occurence')
% set(gca,'FontSize',12)
% 
%     % Crosstalk Ratio promotor strength penalty 
% figure;
% histogram(large_scale_result_file.large_scale_fitting_result_struct.all_crosstalk_ratio_promotor_strength_penalty)
% title('Large Positive Crosstalk Penalty')
% xlabel('Large Positive Crosstalk Penalty')
% ylabel('# Occurence')
% set(gca,'FontSize',12)
% 
%     % Promotor strength sensitivity penalty 
% figure;
% histogram(large_scale_result_file.large_scale_fitting_result_struct.all_promotor_sensitivity_penalty)
% title('Promotor Strength Sensitivity Penalty')
% xlabel('Promotor Strength Sensitivity Penalty')
% ylabel('# Occurence')
% set(gca,'FontSize',12)
% 
%     % mRNA degradataion penalty 
% figure;
% large_scale_result_file.large_scale_fitting_result_struct.all_other_penalty(large_scale_result_file.large_scale_fitting_result_struct.all_other_penalty < 0) = 1e-04; 
% histogram(log10(large_scale_result_file.large_scale_fitting_result_struct.all_other_penalty))
% title('mRNA Degradation Penalty')
% xlabel('log(mRNA Degradation Penalty)')
% ylabel('# Occurence')
% set(gca,'FontSize',12)
%% Implement Sampling using sbiosobol
    
%     % Variations
%         % - Add multiple doses and calculate observables 
%         % - Use stricter bounds for mRNA degradation related parameters 
%         % - Add TX-level crosstalk data 
% 
%         % load fitResult to start with 
% % resultFileName = 'param_est_run_save/20240515_param_est_run0328_objString_15.mat';
% % resultFile = load(resultFileName);
% % fitResult = resultFile.all_fitResults{resultFile.sort_idx(1)}; 
% % problemObject = resultFile.problemObject;
% % modelObj = problemObject.Model; 
% 
%     % Alternatively, use the updated bounds 
% group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
%     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
%     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
%     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
% param_info_path = 'parameters_v2.xlsx';
% problemObject = setProblemObject_v2(group_description,param_info_path,[]); 
% modelObj = problemObject.Model; 
% 
% % Assign estimated values to model parameters
% estimated_param_names = {problemObject.Estimated.Name};
% % estimated_param_values = fitResult.estimated_params; 
% % modified_Mobj = assign_model_parameters(modelObj,estimated_param_names,estimated_param_values); 
% 
% % Add in parameter bounds
% param_bounds = [problemObject.Estimated.Bounds]; 
% param_bounds = reshape(param_bounds,[2,length(estimated_param_names)]);
% param_bounds = param_bounds';
% 
% % Assign input values for sbiosobol
% observables = {'protein sfGFP*'}; 
% 
% % Change dosing information into a vector of doses
% dosing_information = problemObject.Doses;
% 
% %     % Add scenarios to simbiology model 
% % simbiology_scenarios = SimBiology.Scenarios; % Create an empty scenario object 
% % for row_idx = 1:size(dosing_information,1)
% %     add(simbiology_scenarios,'cartesian','dose',dosing_information(row_idx,:)); 
% % end
% 
% % sbiosobol does not support scenarios with dose objects 
%     % Let's try manually assigning the initial values into scenarios 
% simbiology_scenarios = SimBiology.Scenarios; % Create an empty scenario object 
%     % Create target name (which should be the same for all 112 scenarios)
% test_dosing_vector_row = dosing_information(1,:); 
% target_name_list = {test_dosing_vector_row.Target};
% for target_idx = 1:length(target_name_list)
%     target_name = target_name_list{target_idx};
%     target_dosing_vec = dosing_information(:,target_idx); 
%     temp_all_target_dosing = [target_dosing_vec.Amount]; 
%     if ~isempty(temp_all_target_dosing) % If target is involved in dosing
%         all_target_dosing = nan(size(dosing_information,1),1); 
%         for row_idx = 1:size(dosing_information,1)
%             if isempty(dosing_information(row_idx,target_idx).Amount)
%                 all_target_dosing(row_idx) = 0;
%             else
%                 all_target_dosing(row_idx) = dosing_information(row_idx,target_idx).Amount;
%             end
%         end
%         add(simbiology_scenarios,'elementwise',target_name,all_target_dosing); 
%     end
% end
% 
% % Write penalty terms as observables in the model 
% 
% sobolResults = sbiosobol(modelObj,simbiology_scenarios,observables,...
%     'ShowWaitbar',true);

% 
%     % Plot time-course sobol indices 
% timeVec = sobolResults.Time; 
% SobolIndices = sobolResults.SobolIndices; 
% param_names = {SobolIndices.Parameter}; 
% current_mod10 = 0; 
% for param_idx = 1:32
% 
%     % Start a new figure every 10 parameters to keep figure legible 
%     if ~isequal(mod(param_idx,10),0)
%         param_idx_mod10 = mod(param_idx,10);
%     else
%         param_idx_mod10 = 10;
%     end
%     param_idx_divide = ceil(param_idx / 10); 
%     if param_idx_divide > current_mod10
%         figure;
%         current_mod10 = param_idx_divide;
%     end
% 
%     % 
%     param_name = param_names{param_idx}; 
%     FirstOrder = SobolIndices(param_idx).FirstOrder; 
%     TotalOrder = SobolIndices(param_idx).TotalOrder; 
%     if isequal(current_mod10,3)
%         num_row_subplot = 12; 
%     else
%         num_row_subplot = 10; 
%     end
%     subplot(num_row_subplot,2,(param_idx_mod10 - 1) * 2 + 1)
%     plot(timeVec,FirstOrder)
%     hold on 
%     ylabel(strrep(param_name,'_',' '))
%     subplot(num_row_subplot,2,param_idx_mod10 * 2)
%     plot(timeVec,TotalOrder)
% 
% end

%% Implement additional constraints for parameter estimation 
    % Use a result file as a start 

% resultFileName = 'param_est_run_save/20240515_param_est_run0328_objString_15.mat';
% resultFile = load(resultFileName);
% problemObject = resultFile.problemObject;
% 
%     % Apply heuristic bounds on parameter relationships (as Ax <= b) 
% Estimated = problemObject.Estimated;
% param_oi_list_promotor_strength = {'TXTL_PT7_RNAPbound_R','TXTL_PT773_RNAPbound_R','TXTL_PJ23119_RNAPbound_R','TXTL_PJ23105_RNAPbound_R'}; 
% 
%     % 3 inequality constraints: (1) PT7_RNAPbound_R < PT773_RNAPbound_R (2) PJ23119_RNAPbound_R <
%     % PJ23105_RNAPbound_R (3) PT773_RNAPbound_R < PJ23105_RNAPbound_R 
% base_inequality_matrix_promotor_strength = [1 -1 0 0; 0 0 1 -1;0 1 0 -1;]; 
% base_inequality_rhs_promotor_strength = zeros(size(base_inequality_matrix_promotor_strength,1),1); 
% 
% [inequality_A_promotor_strength,inequality_b_promotor_strength] = add_linear_param_constraint(param_oi_list_promotor_strength,base_inequality_matrix_promotor_strength,base_inequality_rhs_promotor_strength,Estimated); 
% 
% param_oi_list_mRNA_deg = {'TXTL_RNAdeg_R_sfGFP','TXTL_RNAdeg_R','TXTL_RNAdeg_kc_sfGFP','TXTL_RNAdeg_kc_kanR','TXTL_RNAdeg_kc'};
% base_inequality_matrix_mRNA_deg = [-1 1 0 0 0;0 0 1 0 -1;0 0 0 1 -1;]; 
% base_inequality_rhs_mRNA_deg = zeros(size(base_inequality_matrix_mRNA_deg,1),1); 
% [inequality_A,inequality_b] = add_linear_param_constraint(param_oi_list_mRNA_deg,base_inequality_matrix_mRNA_deg,...
%     base_inequality_rhs_mRNA_deg,Estimated,inequality_A_promotor_strength,inequality_b_promotor_strength); 

%     % Adding MM fitting 
%         % Use obj string idx = 15 as case study 
% % resultFileName = 'param_est_run_save/20240515_param_est_run0328_objString_15.mat';
% % resultFile = load(resultFileName);
% problemObject = resultFile.problemObject;
% model = problemObject.Model; 
% 
% species_name_list = {model.Species.Name};
% rxn_name_list = {model.Reactions.Reaction};
% remove_idx = strcmp('protein sfGFP*',species_name_list); 
% additional_track_species_list = species_name_list(~remove_idx); 
% simulated_data_name_list = [{'protein sfGFP*'},additional_track_species_list];
% 
% species_name_list = {model.Species.Name}; 
% estimated_params_name_list = {problemObject.Estimated.Name}; 
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% simFunction = resultFile.simFunction; 
% all_fitResults = resultFile.all_fitResults; 
% 
%     % Target reaction rates: mRNA production? 
% for iter = 1:length(all_fitResults)
%     fitResult = all_fitResults{iter}; 
%     estimated_params = fitResult.estimated_params; 
%     [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
%     % rxn_oi = '[CUTP:AGTP:t7RNAP:DNA pT7--utrGFP--sfGFP] -> [term_t7RNAP:DNA pT7--utrGFP--sfGFP] + [RNA utrGFP--sfGFP]'; 
%     rxn_oi = '[CUTP:AGTP:t7RNAP:DNA pT7--utrGFP--sfGFP] -> [term_t7RNAP:DNA pT7--utrGFP--sfGFP] + [RNA utrGFP--sfGFP]'; 
%     rxn_idx = strcmp(rxn_name_list,rxn_oi);
% 
%     all_ss_rxn_rate = nan(length(conc_vec),1); 
%     for conc_idx = 1:length(conc_vec)
% 
%             % Select a baseline case to calculate mRNA production rate 
%         simulated_time_single = simulated_time{conc_idx}; 
%         simulated_data_single = simulated_data{conc_idx};
%         [rxn_name_list,rxn_time_course] = get_reaction_time_course(model,estimated_params_name_list,estimated_params,...
%         simulated_time_single,simulated_data_single,simulated_data_name_list);
% 
%             % Get the reaction rate over time
%         rxn_oi_time_course = rxn_time_course(:,rxn_idx); 
% 
%             % Check for steady state 
% 
%             % Output endpoint rate (assuming ss) 
%         rxn_oi_ss = rxn_oi_time_course(end); 
%         all_ss_rxn_rate(conc_idx) = rxn_oi_ss; 
%     end
% 
%     y = 1 ./ all_ss_rxn_rate; 
%     x = 1 ./ conc_vec; 
%     fitobject = polyfit(x',y,1);
% 
% end
%% Compile results for single objective optimization 

% obj_string_idx_oi = 20; 
% obj_string = all_obj_string{obj_string_idx_oi}; 
%     % 16 - 22 for single objective 
% directory_path = 'param_est_run_save'; 
% directory_files = dir(directory_path); 
% file_names = {directory_files.name};
% 
% missing_iter_list = []; 
% penalty_init = nan(48,1);
% penalty_final = nan(48,1); 
% 
% for iter = 1:48
% 
%     resultFileName_idx = contains(file_names,sprintf('objString_%d_iter_%d',obj_string_idx_oi,iter)); 
%     if any(resultFileName_idx)
%         resultFileName = file_names{resultFileName_idx};    
%         resultFile = load(resultFileName); 
%         penalty_term_idx = strcmp(obj_string,resultFile.penalty_term_label); 
%         penalty_init(iter) = resultFile.penalty_term_init(penalty_term_idx); 
%         penalty_final(iter) = resultFile.penalty_term_final(penalty_term_idx); 
%     else
%         missing_iter_list = [missing_iter_list,iter];
%     end
% 
% end


%% Check result from goal attainment 

% resultFileName = 'param_est_run_save/20240523_param_est_run0954_goal_attainment_run1_iter_1.mat'; 
% resultFile = load(resultFileName); 
% 
% problemObject = resultFile.problemObject;
% dosing_information = create_dosing_info_from_problemObject(problemObject); 
% simFunction = create_simFun_from_problemObject(problemObject); 
% estimated_params = resultFile.estimated_params; 
% [simulated_time,simulated_data] = simFunction(estimated_params,tEnd,dosing_information,tStart:tEnd);
% plot_simulated(simulated_time,simulated_data,'baseline')
% plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')











