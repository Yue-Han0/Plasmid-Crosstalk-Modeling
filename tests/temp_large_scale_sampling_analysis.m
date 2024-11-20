% Perform analysis on parameter sampling results 
clear
clc

currentpath = pwd; 
addpath(genpath(currentpath))
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))

    % Load summary of sampling results 
% result_file = load('test_save_files/20240601_large_scale_parameter_sampling_summary.mat'); 
result_file = load('test_save_files/20240701_sampling_constrained_result_summary.mat'); 
all_sampled_params = result_file.all_sampled_params; 
all_penalty_terms = result_file.all_calculated_penalty_terms; 
% all_dosing_mat = result_file.all_dosing_mat; 
sampled_params_file = load('test_save_files/20240613_param_init_cond_sampling_constrained.mat'); 
all_dosing_information = sampled_params_file.all_dosing_information; 
% sampled_params_file = load('test_save_files/20240701_sam')


    % Load an example result file to get parameter names 
group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
    'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
    'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
    'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
param_info_path = 'parameters_v2.xlsx';
problemObject = setProblemObject_v2(group_description,param_info_path,[]);
param_names = {problemObject.Estimated.Name}; 

% Preprocess to remove missing values if necessary? 
invalid_row_idx_list = any(isnan(all_penalty_terms),2); 
valid_all_penalty_terms = all_penalty_terms(~invalid_row_idx_list,:); 
valid_all_sampled_params = all_sampled_params(~invalid_row_idx_list,:); 
% 
%     % Let's do a sanity check that the sampled parameters in the result
%     % analysis file matches with the selected parameters in sampled params
%     % file 
% sampled_params_file_nominal_sampled_params = sampled_params_file.nominal_sampled_params;
% sampled_params_file_nominal_sampled_params_selected = sampled_params_file_nominal_sampled_params(sampled_params_file.satisfied_idx_list,:);
% result_file_sampled_params = all_sampled_params(1:100:end,:); 
% not_matched_idx_list = find(~all(abs(sampled_params_file_nominal_sampled_params_selected - result_file_sampled_params(1:9578,:)) < 1e-04,2));
% unmatched_sampled_params_file_params = sampled_params_file_nominal_sampled_params_selected(not_matched_idx_list,:); 
% unmatched_result_file_params = result_file_sampled_params(not_matched_idx_list,:); 


%% 1.1 Plot distribution of penalty terms 
% for penalty_idx = 1:size(all_penalty_terms,2)
%     penalty_term_label = result_file.penalty_labels{penalty_idx};  
%     penalty_term = valid_all_penalty_terms(:,penalty_idx); 
%     % Remove outliers 
%     % filtered_penalty_term = rmoutliers(penalty_term); 
%     filtered_penalty_term = penalty_term; 
% 
%     if mean(penalty_term) > 10 && all(filtered_penalty_term > 0)
%         penalty_term_for_plot = log10(filtered_penalty_term);
%         xlabel_name = sprintf('log10(%s)',strrep(penalty_term_label,'_',' ')); 
%     else
%         penalty_term_for_plot = filtered_penalty_term;
%         xlabel_name = strrep(penalty_term_label,'_',' '); 
%     end
%     figure; 
%     histogram(penalty_term_for_plot)
%     % histogram(penalty_term_for_plot(penalty_term_for_plot < 0))
%     xlabel(xlabel_name)
%     ylabel('# Occurence')
%     title(strrep(penalty_term_label,'_',' ')); 
%     set(gca,'FontSize',12);
%     saveas(gcf,sprintf('plots/20240723_sampling_result_summary0701_%s_distribution.png',penalty_term_label)); 
% 
% end

%% 1.2 parameter distribution for satisfied vs. unsatisfeid for each penalty
% satisfied_penalty_cut_off_list = [nan,nan,0,0,-1,1e+04]; 
% for penalty_idx = 3:size(all_penalty_terms,2)
% 
%     % Get penalty term labels for plot titles & cut-off values 
%     penalty_term_label = result_file.penalty_labels{penalty_idx}; 
%     satisfied_penalty_cut_off = satisfied_penalty_cut_off_list(penalty_idx); 
% 
%     figure; 
%     for param_idx = 1:size(valid_all_sampled_params,2)
%         param_name = param_names{param_idx}; 
%         subplot(ceil(sqrt(size(valid_all_sampled_params,2))),ceil(sqrt(size(valid_all_sampled_params,2))),param_idx)
% 
%         % Separate sampled parameters into satisfied vs. unsatisfied group 
%         satisfied_sampled_penalty_term_idx = valid_all_penalty_terms(:,penalty_idx) < satisfied_penalty_cut_off;
%         unsatisfied_sampled_penalty_term_idx = valid_all_penalty_terms(:,penalty_idx) > satisfied_penalty_cut_off;
% 
%         satisfied_sampled_params = valid_all_sampled_params(satisfied_sampled_penalty_term_idx,param_idx);
%         unsatisfied_sampled_params = valid_all_sampled_params(unsatisfied_sampled_penalty_term_idx,param_idx);
% 
%         % Plot out distribution in log space
%         histogram(log10(satisfied_sampled_params))
%         hold on 
%         histogram(log10(unsatisfied_sampled_params))
%         title(strrep(param_name,'_',' '))
% 
%     end
%     legend('Satisfied','Unsatified')
%     sgtitle(strrep(penalty_term_label,'_',' '))
% end 


%% 2.1 Correlation-based sensitivity analysis. using PCC & PRCC - with parameters and initial conditions
% separately 
    % Note that these are recommended when the relationship is monotonic,
    % which may or may not be true here 

% PRCC = partialcorri(valid_all_penalty_terms,valid_all_sampled_params); 
% figure; 
% h1 = heatmap(strrep(param_names,'_',' '),strrep(result_file.penalty_labels,'_',' '),PRCC,'Colormap',jet); 
% xlabel('Parameter Names')
% ylabel('Penalty Terms') 
% 
% % Let's try using PCC & PRCC - with both parameters and initial conditions 
% mega_dosing_mat = repmat(sampled_params_file.nominal_sampled_empty_init_cond,10000,1); 
% valid_mega_dosing_mat = mega_dosing_mat(~invalid_row_idx_list,:);
% augmented_param_initCond_mat = [valid_all_sampled_params,valid_mega_dosing_mat];
% augmented_param_initCond_labels = [strrep(param_names,'_',' '),{'Empty Plasmid Concentration'}]; % Param names and initial cond names 
% augmented_PCC = partialcorri(valid_all_penalty_terms,augmented_param_initCond_mat); 
% rank_valid_all_penalty_terms = tiedrank(valid_all_penalty_terms);
% rank_augmented_param_initCond_mat = tiedrank(augmented_param_initCond_mat); 
% augmented_PRCC = partialcorri(rank_valid_all_penalty_terms,rank_augmented_param_initCond_mat); 
% 
% figure; 
% h2 = heatmap(augmented_param_initCond_labels,strrep(result_file.penalty_labels,'_',' '),augmented_PCC,'Colormap',jet); 
% xlabel('Parameter & Initial Condition Names')
% ylabel('Penalty Terms') 
% title('Partial Correlation Coefficient')
% 
% figure; 
% h3 = heatmap(augmented_param_initCond_labels,strrep(result_file.penalty_labels,'_',' '),augmented_PRCC,'Colormap',jet); 
% xlabel('Parameter & Initial Condition Names')
% ylabel('Penalty Terms') 
% title('Partial Rank Correlation Coefficient')
% 
% % What about calculating just correlation? 
% corr_matrix = corr(valid_all_penalty_terms,augmented_param_initCond_mat); 
% figure;
% h4 = heatmap(augmented_param_initCond_labels,strrep(result_file.penalty_labels,'_',' '),corr_matrix,'Colormap',jet); 
% xlabel('Parameter & Initial Condition Names')
% ylabel('Penalty Terms')


%% 3. Is there a sample that satisfies all phenotype criteria? 
% 
%     % We don't really care about deviations 
% satisfied_idx_list = valid_all_penalty_terms(:,3) < 0 & valid_all_penalty_terms(:,4) < 0 & valid_all_penalty_terms(:,5) < -1 & valid_all_penalty_terms(:,6) < 1e+04; 
%     % Check the sampled parameter values and initial condition 
% satisfied_sampled_params = valid_all_sampled_params(satisfied_idx_list,:); 
% satisfied_init_cond = valid_mega_dosing_mat(satisfied_idx_list,:); 
% satisfied_penalty_terms = valid_all_penalty_terms(satisfied_idx_list,:); 
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
%     histogram(log10(selected_satisfied_sampled_params))
%     title(strrep(param_name,'_',' '));
%     xlabel('log_{10}parameter')
% end
% sgtitle('Parameter Distribution Satisfying Qualitative Phenotypes');

    % What does the deviation look like? 
% satisfied_penalty_terms = valid_all_penalty_terms(satisfied_idx_list,:);

%     % Plot out 10 of these 
%         % Create SimFunction first
%     % Calculate time-course data using the current set of parameters 
% dosing_target_list = sampled_params_file.target_name_list;
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
%     sampled_params_oi = satisfied_sampled_params(selected_idx,:); 
%     sampled_penalties_oi = satisfied_penalty_terms(selected_idx,:); 
%     sampled_init_cond_oi = satisfied_init_cond(selected_idx,:);
%     init_cond_oi_idx = find(sampled_params_file.nominal_sampled_empty_init_cond == sampled_init_cond_oi); 
%     dosing_info_oi = all_dosing_information{init_cond_oi_idx}; 
% 
%     [simulated_time,simulated_data] = simFunction(sampled_params_oi,tEnd,dosing_info_oi,tStart:tEnd);
%     plot_simulated(simulated_time,simulated_data,'baseline')
%     plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')
% 
% end
% 
%     % Filter out cases where large positive crosstalk ratio is satisfied 
% large_positive_crosstalk_ratio_penalty = all_penalty_terms(:,5); 
% satisfied_idx_list = find(large_positive_crosstalk_ratio_penalty < -1); % note that the threshold is arbitrarily set to -1 for now  
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

%% 4. Correlation with combinations of parameters 
    % Calculate correlation of penalty terms with pairwise parameters
    % multiplied or added (i.e. mimicing kcat * E0 or (k-1 + kcat)/k1

    % Preassign 
param_combo_name_list = cell(930,1); 
param_combo_value_list = nan(size(valid_all_sampled_params,1),930);  

combo_idx = 1; 
for param_idx_1 = 1:length(param_names) - 1

    % Get parameter A
    param_name_1 = param_names{param_idx_1}; 
    param_val_1 = valid_all_sampled_params(:,param_idx_1); 
    for param_idx_2 = param_idx_1 + 1:length(param_names)

        % Get parameter B 
        param_name_2 = param_names{param_idx_2}; 
        param_val_2 = valid_all_sampled_params(:,param_idx_2); 

        % Transform
        param_combo_add_name = sprintf('%s+%s',param_name_1,param_name_2); 
        param_combo_multi_name = sprintf('%s*%s',param_name_1,param_name_2); 
        param_combo_add = param_val_1 + param_val_2; 
        param_combo_multi = param_val_1 .* param_val_2; 

        % Store 
        param_combo_name_list{combo_idx} = param_combo_add_name;
        param_combo_name_list{combo_idx + 1} = param_combo_multi_name;
        param_combo_value_list(:,combo_idx) = param_combo_add; 
        param_combo_value_list(:,combo_idx + 1) = param_combo_multi; 

        % Update index 
        combo_idx = combo_idx + 2; 

    end
end

    % Remove duplicate columns 
[final_param_combo_name_list,param_combo_unique_idx] = unique(param_combo_name_list);
final_param_combo_value_list = param_combo_value_list(:,param_combo_unique_idx); 

%     % 
% [PRCC_param_combo,param_combo_pval] = partialcorri(valid_all_penalty_terms,final_param_combo_value_list);

    % Calculate partial rank correlation 
rank_valid_all_penalty_terms = tiedrank(valid_all_penalty_terms);
rank_param_combo_mat = tiedrank(final_param_combo_value_list); 
[PRCC_param_combo,param_combo_pval] = partialcorri(rank_valid_all_penalty_terms,rank_param_combo_mat); 

save('test_save_files/20240729_param_combo_partial_rank_correlation.mat'); 
%     % Apply Benjamini-Hochberg procedure to correct these p-values 
% all_FDR = nan(size(param_combo_pval)); 
% for penalty_idx = 1:size(param_combo_pval,1)
%     FDR = mafdr(param_combo_pval(penalty_idx,:));
%     all_FDR(penalty_idx,:) = FDR; 
% end
% 
% 
% 
%     % Use a 5% false discovery rate as the threshold 
% keep_col_idx_list = any(all_FDR < 0.05,1);
% significant_PRCC = PRCC_param_combo(:,keep_col_idx_list); 
% significant_param_combo_names = final_param_combo_name_list(keep_col_idx_list); 
% % significant_PRCC(all_FDR >= 0.05) = 0;
% heatmap_ylabels = strrep(result_file.penalty_labels,'_',' ');
% heatmap_xlabels = strrep(significant_param_combo_names,'_',' '); 
% figure; 
% heatmap(heatmap_xlabels,heatmap_ylabels,significant_PRCC); 
% title('Statistically Significantly Correlated Parameter Combinations')
% 
%     % Also just check correlation coefficients whose absolute values are larger than 0.15
% keep_col_idx_list_2 = any(abs(PRCC_param_combo) > 0.15); 
% significant_PRCC_2 = PRCC_param_combo(:,keep_col_idx_list_2); 
% significant_param_combo_names_2 = final_param_combo_name_list(keep_col_idx_list_2); 
% % significant_PRCC(all_FDR >= 0.05) = 0;
% heatmap_ylabels_2 = strrep(result_file.penalty_labels,'_',' ');
% heatmap_xlabels_2 = strrep(significant_param_combo_names_2,'_',' '); 
% figure; 
% heatmap(heatmap_xlabels_2,heatmap_ylabels_2(1:end-1),significant_PRCC_2(1:end-1,:)); 
% title('Correlated Parameter Combinations (abs value > 0.15')


%% 5. Is the current number of samples enough? 
% N = 300000; 
% num_draft = 10; 
% 
% drafted_bs_sampled_params = cell(num_draft,1);
% drafted_bs_penalty_terms = cell(num_draft,1); 
% 
%     % Draft parameters and penalty terms 
% for draw_idx = 1:num_draft
%         % Draw N samples 
%     lottery = randperm(size(valid_all_penalty_terms,1)); 
%     select_idx_list = lottery(1:N); 
% 
%         % Get bootstrapped sampled params and penalty terms 
%     bs_sampled_params = valid_all_sampled_params(select_idx_list',:);
%     bs_penalty_terms = valid_all_penalty_terms(select_idx_list,:); 
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
%     sgtitle(strrep(result_file.penalty_labels{penalty_idx},'_',' '))
% 
% end



 
