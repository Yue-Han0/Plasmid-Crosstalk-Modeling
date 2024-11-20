clear   
clc

addpath(genpath('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling'));
rmpath('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/benchmark_model')

%% Test crosstalk prediction with default parameters 
group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
    'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
    'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
param_info_path = '/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/benchmark_model/parameters_ori.xlsx';
problemObject = setProblemObject_updated(group_description,param_info_path);
Mobj = problemObject.Model; 

% Plot crosstalk with untrained model (confirm that crosstalk is there) 
    % Change initial concentrations of RNAP/RNase
    RNAP = sbioselect(Mobj,'Type','Parameter','Name','RNAP_0');
    t7RNAP = sbioselect(Mobj,'Type','Parameter','Name','t7RNAP_0');
    RNase = sbioselect(Mobj,'Type','Parameter','Name','RNase_0');
    Ribo = sbioselect(Mobj,'Type','Parameter','Name','Ribo_0');
    set(t7RNAP,'Value',7.5);
    set(RNase,'Value',20); 
    % set(Ribo,'Value',200);

    pT7_RNAPbound_R = sbioselect(Mobj,'Type','Parameter','Name','TXTL_PT7_RNAPbound_R'); %0.1
    pT773_RNAPbound_R = sbioselect(Mobj,'Type','Parameter','Name','TXTL_PT773_RNAPbound_R');
    pJ23119_RNAPbound_R = sbioselect(Mobj,'Type','Parameter','Name','TXTL_PJ23119_RNAPbound_R');
    pkanR_RNAPbound_R = sbioselect(Mobj,'Type','Parameter','Name','TXTL_PkanR_RNAPbound_R');
        % Larger value for weaker promotors 
    set(pT773_RNAPbound_R,'Value',100);
    set(pJ23119_RNAPbound_R,'Value',10);
    set(pkanR_RNAPbound_R,'Value',1);
    
    [test_simData] = txtl_runsim(Mobj,14 * 60 * 60);
    % Change plasmid concentration and plot out crosstalk 
    all_promotor_description = {'T7_strong','T7_weak','sigma70_strong'}; 
    all_empty_description = {'no_empty','empty','empty_T7','empty_sigma70'};
    all_plasmid_concentration = [0.5,1,2.5,5,10,15,30]; 
    all_colors = {'g',[0.5,0.5,0.5],'r','b'};
    all_markers = {'o','+','x','square'}; 
    figure;
    for conc_idx = 1:length(all_plasmid_concentration) % Assume 7 sets of data are always used 
        plasmid_conc = all_plasmid_concentration(conc_idx); 
        for prom_idx = 1:length(all_promotor_description)

            plot_idx = (conc_idx - 1) * length(all_promotor_description) + prom_idx;
            subplot(length(all_plasmid_concentration),length(all_promotor_description),plot_idx)
            
            promotor_description = all_promotor_description{prom_idx};
            if isequal(conc_idx,1)
                title(strrep(promotor_description,'_',' '))
                hold on 
            end
            if isequal(conc_idx,7)
                xlabel('time(s)')
                hold on 
            end
            if isequal(prom_idx,1)
                ylabel(sprintf('%.1f nM',plasmid_conc))
                hold on 
            end
            for empty_idx = 1:length(all_empty_description)

                temp_Mobj = copyobj(Mobj); 
                empty_description = all_empty_description{empty_idx}; 
                desired_data_description = strcat(promotor_description,'_',empty_description,'_',num2str(plasmid_conc)); 
                
                % Change initial condition for plasmids 
                    % Reporter 
                if contains(desired_data_description,'T7_strong')
                    promotor_name = 'T7';
                elseif contains(desired_data_description,'T7_weak')
                    promotor_name = 'T773'; 
                elseif contains(desired_data_description,'sigma70_strong')
                    promotor_name = 'J23119';
                end
                reporter_dna = sbioselect(temp_Mobj,'Type','Species','Name',['DNA p' promotor_name '--utrbroc--no_protein']);
                set(reporter_dna,'Value',plasmid_conc)
                    % empty 
                if contains(desired_data_description,'no_empty')
                    empty_name = '';
                elseif contains(desired_data_description,'empty_T7')
                    empty_name = 'T7';
                elseif contains(desired_data_description,'empty_sigma70')
                    empty_name = 'J23119'; 
                elseif contains(desired_data_description,'empty')
                    empty_name = 'empty';
                end
                if ~isempty(empty_name)
                    empty_conc = 15;
                end
                empty_dna = sbioselect(temp_Mobj,'Type','Species','Name',['DNA p' empty_name '--utrempty--no_protein']);
                if ~isempty(empty_dna)
                    set(empty_dna,'Value',empty_conc); 
                end
                
                    % kanR 
                kan_dna = sbioselect(temp_Mobj,'Type','Species','Name','DNA pkanR--utrkanR--kanR');
                if exist('empty_conc','var')
                    set(kan_dna,'Value',plasmid_conc + empty_conc)
                else
                    set(kan_dna,'Value',plasmid_conc)
                end

                [simData] = txtl_runsim(temp_Mobj,14 * 60 * 60); 
                predicted_mRNA_conc_idx = strcmp(simData.DataNames,'RNA utrbroc--no_protein');
                predicted_mRNA_conc = simData.Data(:,predicted_mRNA_conc_idx); 
                plot(simData.Time,predicted_mRNA_conc,'LineWidth',1.5,'Marker',all_markers{empty_idx},'Color',all_colors{empty_idx})
                hold on

            end
        end
    end

line1 = plot(nan,nan,'Color','g','LineStyle','--','DisplayName','no empty simulated');
line2 = plot(nan,nan,'Color',[0.5,0.5,0.5],'LineStyle','--','DisplayName','empty simulated');
line3 = plot(nan,nan,'Color','r','LineStyle','--','DisplayName','empty T7 simulated');
line4 = plot(nan,nan,'Color','b','LineStyle','--','DisplayName','empty sigma70 simulated');
legend([line1,line2,line3,line4])




%% Test txtl_sim_predict
% fit_result_path = 'param_est_run_save/20230413_param_est_run1309_modify_elongation_rate.mat';
% fit_result_path = 'param_est_run_save/20230413_param_est_run1915_.mat'; 
% group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% txtl_sim_predict(fit_result_path,group_description)

%% Calculate CI for fitted model
% test_fitResultFile = load('param_est_run_save/20230321_sigma70_param_est_run1129.mat'); 
% test_fitResult_all = test_fitResultFile.all_fitResults;
% test_fitResult = test_fitResult_all{1}; 
% [simulated_data,param_estimate] = predict(test_fitResult); 
% [fitted_data,param_estimate2] = fitted(test_fitResult); 
% 
% numRuns = 48; 
% % jacobian = test_fitResult.J; 
% % allZero_col_idx_list = all(jacobian==0); 
% % corresponding_param_list = test_fitResultFile.problemObject.Estimated(allZero_col_idx_list);
% % revised_jacobian = jacobian(:,~allZero_col_idx_list); 
% % jacobian_sq = revised_jacobian' * revised_jacobian; 
% % revised_covb = inv(revised_jacobian' * revised_jacobian) * test_fitResult.MSE; 
% % revised_standardError = sqrt(diag(revised_covb)); % In log space 
% % % update test_fitResult
% %     % Beta,
% %     % ParameterEstimates,J,COVB,CovarianceMatrix,EstimatedParameterNamees 
% % skip_idx_list = find(allZero_col_idx_list);
% % revised_test_fitResult = test_fitResult; 
% %     % update StandardError in Beta & ParameterEstimates
% % revised_error_idx = 1;
% % for param_idx = 1:height(test_fitResult.Beta)
% %     if ~ismember(param_idx,skip_idx_list)
% %         table_idx_name = sprintf('beta(%d)',param_idx);
% %         revised_test_fitResult.Beta.StandardError(table_idx_name) = revised_standardError(revised_error_idx); 
% %         revised_error_idx = revised_error_idx + 1; 
% %         % nominal standard error = exp(log mean + log standard error) -
% %         % exp(log mean) 
% %         nominal_standard_error = exp(test_fitResult.Beta.Estimate(table_idx_name) +...
% %             revised_standardError(revised_error_idx)) - exp(test_fitResult.Beta.Estimate(table_idx_name)); 
% %         revised_test_fitResult.ParameterEstimates.StandardError(param_idx) = nominal_standard_error; 
% %     end
% % 
% % end
% % Get SSE and estiamted parameter summary 
% SSE_summary = nan(numRuns,1); 
% log_likelihood_summary = nan(numRuns,1); 
% 
% for iter = 1:numRuns
%     fitResult = all_fitResults{iter,1};  
%     SSE_summary(iter,1) = fitResult.SSE;
%     log_likelihood_summary(iter,1) = fitResult.LogLikelihood; 
%     est_param_table = fitResult.ParameterEstimates; 
%     if ~exist('estimated_param_summary','var')
%         estimated_param_summary = est_param_table;
%         estimated_param_summary = removevars(estimated_param_summary,'StandardError');
%         estimated_param_summary = renamevars(estimated_param_summary,'Estimate','Estimate_1');
%     else
%         new_col_name = sprintf('Estimate_%d',iter); 
%         estimated_param_summary.(new_col_name) = est_param_table.Estimate; 
%     end
% 
% end
% paramci = sbioparameterci(test_fitResult); 
% predci = sbiopredictionci(opt_fitResult); 
    % Where do the nan values in residuals come from? from the simulated
    % data? from experimental data? 
% Sanity check #1: Does the #timepoints in simulated data add up to number
% of residuals? 
% total_num_timepoints = 0; 
% for i = 1:length(fitted_data)
%     fitted_data_obj = fitted_data(i); 
%     rna_idx = strcmp(fitted_data_obj.DataNames,'RNA utrbroc--no_protein'); 
%     rna_time_course = fitted_data_obj.Data(:,rna_idx); 
%     num_timepoint = length(rna_time_course); 
%     total_num_timepoints = total_num_timepoints + num_timepoint; 
% end
% It does not. The #timepoints corresponds to the experimental data 

%% Plot baseline fitting results with the new data storage format & confidence intervals 
% test_fitResultFile = load('param_est_run_save/20230307_sigma70_param_est_run1251.mat'); 
% test_fitResult_all = test_fitResultFile.all_fitResults;
% 
% % Get SSE summary 
% numRuns = 49; 
% SSE_summary = nan(numRuns,1); 
% for iter = 1:numRuns
%     fitResult = test_fitResult_all{iter + 1,1};  
%     SSE_summary(iter,1) = fitResult.SSE;
% end
% [~,opt_idx] = min(SSE_summary); 
% opt_fitResult = test_fitResult_all(opt_idx); 
% [fitted_data,param_estimate2] = fitted(opt_fitResult{1,1});
% 
% group_description = {'sigma70_no_empty_lowConc','sigma70_no_empty_highConc'};
% param_info_path = '/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/benchmark_model/parameters_test1.xlsx';
% 
% problemObject = setProblemObject(group_description,param_info_path); 
% experimental_data = problemObject.Data; 
% 
% figure;
% 
% % Plot experimental data 
% exp_timeVec_all = experimental_data.Time;
% exp_mRNA_conc_all = experimental_data.mRNA_concentration; 
% start_idx_list = find(exp_timeVec_all < 30); 
% for idx = 1:length(start_idx_list)
%     start_idx = start_idx_list(idx);
%     if isequal(idx,length(start_idx_list))
%         end_idx = length(exp_mRNA_conc_all);
%     else
%         end_idx = start_idx_list(idx + 1); 
%     end
%     plot(exp_timeVec_all(start_idx:end_idx),exp_mRNA_conc_all(start_idx:end_idx),'LineStyle','--','LineWidth',1.5)
%     hold on 
% end
% 
% % Plot fitted data
% for i = 1:length(fitted_data)
%     sim_data = fitted_data(i).Data;
%     rna_idx = strcmp(fitted_data(i).DataNames,'RNA utrbroc--no_protein'); 
%     rna_time_course = sim_data(:,rna_idx); 
%     timeVec = fitted_data(i).Time;
%     plot(timeVec,rna_time_course,'LineStyle','-','LineWidth',1.5);
%     hold on 
% end
% legend('0.5nM experimental','1nM experimental',...
%     '2.5nM experimental','5nM experimental','10nM experimental',...
%     '15nM experimental','30nM experimental','0.5nM simulated','1nM simulated',...
%     '2.5nM simulated','5nM simulated','10nM simulated',...
%     '15nM simulated','30nM simulated')
% xlabel('Time (s)')
% ylabel('mRNA concentration')
%% Calculate SSE for experimental data CI boundary 
% 
% load data_structures/simbio_data_table_updated_AP.mat
% % group_description = {'sigma70_strong_no_empty','sigma70_strong_empty',...
% %     'sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% group_description = {'sigma70_strong_no_empty'};
% [group_number,selected_table] = get_relevant_data(all_data_table,variable_name_stem_list,group_description);
% selected_table = remove_neg_val(selected_table); 
% 
% gain = 100; 
% mRNA_concentration = fluo_mRNA_conversion(selected_table.fluorescence,[],gain);
% selected_table = addvars(selected_table,mRNA_concentration,'NewVariableNames','mRNA_concentration'); 
% mRNA_CI_lb = fluo_mRNA_conversion(selected_table.CI_lb,[],gain);
% selected_table = addvars(selected_table,mRNA_CI_lb,'NewVariableNames','mRNA_CI_lb'); 
% mRNA_CI_ub = fluo_mRNA_conversion(selected_table.CI_ub,[],gain);
% selected_table = addvars(selected_table,mRNA_CI_ub,'NewVariableNames','mRNA_CI_ub'); 
%  
% noise_SSE = sum((selected_table.mRNA_CI_lb - selected_table.mRNA_concentration).^2,'omitnan');
% 
% 
% function table_for_fit = remove_neg_val(table_for_fit)
%     for i = 1:length(table_for_fit.fluorescence)
%         if table_for_fit.fluorescence(i) < 0
%             table_for_fit.fluorescence(i) = nan;
%         end
% 
%     end
% 
% end
% function [group_number,selected_table] = get_relevant_data(data_table,data_description,group_description)
% %     group_description - a cell array with description of data to be
% %     selected
% %     group_number - a vector of the group number for selected data
% %     selected_table - a section of data table requested 
% 
%     group_number = []; 
%     for group_idx = 1:length(group_description)
%         group_description_s = group_description{group_idx};
%         group_number = [group_number;find(strcmp(group_description_s,data_description))]; 
%     end
%     selected_idx = []; 
%     for group_num = group_number'
%         selected_idx = [selected_idx;find(data_table.Group==group_num)];
%     end
%     selected_table = data_table(selected_idx,:); 
% 
% end