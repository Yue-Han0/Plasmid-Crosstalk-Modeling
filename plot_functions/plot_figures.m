%
clear
clc

currentpath = pwd; 
addpath(genpath(currentpath))
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))

%% Load result files 

    % Define simulation time span 
    tStart = 0; 
    tEnd = 21600; 

    % Define initial plasmid concentrations used in experiments 
    conc_vec = [0.5,1,2.5,5,10,15,30]; 

    % Define promotor name list and plasmid composition 
    promotor_name_list = {'T7_strong ','T7_weak ','sigma70_strong ','sigma70_weak '};
    composition_name_list = {'empty','empty T7','empty \sigma70'}; 

    % Define names for the 9 criteria
    criteria_name_list = {'Baseline Deviation','Crosstalk Ratio Deviation','Baseline Penalty Deviation','T7 Strong Crosstalk',...
        'T7 weak Crosstalk','sigma70 strong Crosstalk','sigma70 weak Crosstalk','Larger Positive Crosstalk Penalty',...
        'Residual mRNA Penalty'};

    % Load sampled parameters file 
    init_param_file = load('result_files/20240819_lhs_sampled_params.mat'); 
    all_satisfied_sampled_params = init_param_file.satisfied_sampled_params;

    % Load result file 
    result_summary_file = load('result_files/202409_high_res_penalty_sampling_result_summary.mat'); 

    % Load sampled parameter set of interest 
    model_oi = load('result_files/20241023_model_oi_info.mat'); 
    simFunction = model_oi.simFunction;
    problemObject = model_oi.problemObject;
    dosing_information = model_oi.dosing_information; 
    sampled_params_oi = model_oi.sampled_params_oi; 

    % Load penalty cutoff 
    high_res_penalty_goals = load('result_files/high_resolution_penalty_goals.mat'); 
    penalty_cutoff = high_res_penalty_goals.penalty_goals; 

    % Get useful variables 
    species_name_list = simFunction.Observables.Name; 
    parameters_name_list = {problemObject.Estimated.Name};

    % Load an example result file 
    sample_result_file = load('result_files/20240821_param_sampling_constrained_run1.mat');
    high_res_penalty_term_labels = sample_result_file.penalty_term_labels; 
    high_res_penalty_term_length = sample_result_file.penalty_term_length; 

    % Define colors 
    cmap = hsv(256);
    colorIndices = round(linspace(1, size(cmap, 1) - 64,length(conc_vec)));
    sampledColors = cmap(colorIndices, :);

%% Figure 2 

    % % A. Baseline expression & B. Crosstalk ratio
    [simulated_time,simulated_data] = simFunction(sampled_params_oi,tEnd,dosing_information,tStart:tEnd);
    % plot_simulated(simulated_time,simulated_data,'baseline')
    % plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')

%% Figure 3 & S12

    % % Select datasets
    % start_idx = 56; 
    %     % Plasmid conc = 2.5nM 
    % sigma70_strong_lowConc1_no_empty_simTime = simulated_time{59}; 
    % sigma70_strong_lowConc1_empty_simTime = simulated_time{66}; 
    % sigma70_strong_lowConc1_empty_T7_simTime = simulated_time{73}; 
    % sigma70_strong_lowConc1_empty_sigma70_simTime = simulated_time{80}; 
    % sigma70_strong_lowConc1_no_empty_simData = simulated_data{59}; 
    % sigma70_strong_lowConc1_empty_simData = simulated_data{66}; 
    % sigma70_strong_lowConc1_empty_T7_simData = simulated_data{73}; 
    % sigma70_strong_lowConc1_empty_sigma70_simData = simulated_data{80}; 
    % % Plasmid conc = 0.5nM 
    % sigma70_strong_lowConc2_no_empty_simTime = simulated_time{57}; 
    % sigma70_strong_lowConc2_empty_simTime = simulated_time{64}; 
    % sigma70_strong_lowConc2_empty_T7_simTime = simulated_time{71}; 
    % sigma70_strong_lowConc2_empty_sigma70_simTime = simulated_time{78}; 
    % sigma70_strong_lowConc2_no_empty_simData = simulated_data{57}; 
    % sigma70_strong_lowConc2_empty_simData = simulated_data{64}; 
    % sigma70_strong_lowConc2_empty_T7_simData = simulated_data{71}; 
    % sigma70_strong_lowConc2_empty_sigma70_simData = simulated_data{78}; 
    %     % Plasmid conc = 10 nM 
    % sigma70_strong_highConc1_no_empty_simTime = simulated_time{61}; 
    % sigma70_strong_highConc1_empty_simTime = simulated_time{68}; 
    % sigma70_strong_highConc1_empty_T7_simTime = simulated_time{75}; 
    % sigma70_strong_highConc1_empty_sigma70_simTime = simulated_time{82}; 
    % sigma70_strong_highConc1_no_empty_simData = simulated_data{61}; 
    % sigma70_strong_highConc1_empty_simData = simulated_data{68}; 
    % sigma70_strong_highConc1_empty_T7_simData = simulated_data{75}; 
    % sigma70_strong_highConc1_empty_sigma70_simData = simulated_data{82}; 

    % % A. Toxin mechanism diagnosis 
    % figure; 
    % colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    %     [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
    %     [0.6350 0.0780 0.1840]}; 
    % species_oi_list_toxin_mech = {'toxin','Ribo'}; % Check related species concentration over time 
    % parameter_oi_name = 'toxin_threshold'; % Get the threhold value f
    % parameter_oi_idx = find(strcmp(parameter_oi_name,parameters_name_list)); 
    % for species_idx = 1:length(species_oi_list_toxin_mech)
    %     species_oi_name = species_oi_list_toxin_mech{species_idx};
    %     species_oi_idx_list = find(contains(species_name_list,species_oi_name)); 
    %     subplot(ceil(sqrt(length(species_oi_list_toxin_mech))),ceil(sqrt(length(species_oi_list_toxin_mech))),species_idx)
    %     for conc_idx = 1:2:7
    %             % plot no empty
    %         simulated_time_single_oi = simulated_time{start_idx + conc_idx};
    %         simulated_data_single_oi = simulated_data{start_idx + conc_idx};
    %         plot(simulated_time_single_oi ./ 3600,...
    %             sum(simulated_data_single_oi(:,species_oi_idx_list),2),'LineWidth',1.5,'Color',colors{conc_idx})
    %         hold on 
    %             % plot empty with dashed lines 
    %         simulated_time_single_oi_empty = simulated_time{start_idx + length(conc_vec) + conc_idx};
    %         simulated_data_single_oi_empty = simulated_data{start_idx + length(conc_vec) + conc_idx};
    %         plot(simulated_time_single_oi_empty ./ 3600,...
    %             sum(simulated_data_single_oi_empty(:,species_oi_idx_list),2),'LineWidth',1.5,'LineStyle','--','Color',colors{conc_idx})
    %     end
    %     if isequal(species_idx,1)
    %         yline(sampled_params_oi(parameter_oi_idx))
    %     end
    %     xlabel('Time(hr)')
    %     ylabel('Concentration (nM)')
    %     if isequal(species_idx,length(species_oi_list_toxin_mech))
    %         legend('0.5nM baseline','0.5nM +empty','2.5nM baseline','2.5nM +empty',...
    %             '10nM baseline','10nM +empty','30nM baseline','30nM +empty')
    %     end
    %     title(species_oi_name)
    %     set(gca,'FontSize',14)
    % end
    
    % % B. Total concentrations 
    % figure; 
    % 
    % species_oi_list_total_conc = {'RNA utrGFP','RNA utrkanR','RNA utrempty','protein sfGFP*','protein kanR'};
    % species_oi_plot_name_list = {'sfGFP mRNA','kanR mRNA','empty mRNA','sfGFP protein','kanR protein'};
    % for species_idx = 1:length(species_oi_list_total_conc)
    %     species_oi_name = species_oi_list_total_conc{species_idx};
    %     species_oi_idx_list = find(contains(species_name_list,species_oi_name) & ~contains(species_name_list,'RNase'));
    % 
    %     subplot(ceil(sqrt(length(species_oi_list_total_conc))),ceil(sqrt(length(species_oi_list_total_conc))),species_idx)
    %     % Plot low concentration 
    %     plot(sigma70_strong_lowConc2_no_empty_simTime ./ 3600,sum(sigma70_strong_lowConc2_no_empty_simData(:,species_oi_idx_list),2),'LineWidth',1.5,'Color','g')
    %     hold on 
    %     plot(sigma70_strong_lowConc2_empty_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_simData(:,species_oi_idx_list),2),...
    %         'LineWidth',3,'Color',[0.5,0.5,0.5])
    %     plot(sigma70_strong_lowConc2_empty_T7_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_T7_simData(:,species_oi_idx_list),2),...
    %         'LineWidth',2.25,'Color','r')
    %     plot(sigma70_strong_lowConc2_empty_sigma70_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_sigma70_simData(:,species_oi_idx_list),2),...
    %         'LineWidth',1.5,'Color','b')
    % 
    %     % % Plot high concentration 
    %     % plot(sigma70_strong_highConc1_no_empty_simTime ./ 3600,sum(sigma70_strong_highConc1_no_empty_simData(:,species_oi_idx_list),2),'LineWidth',1.5,'Color','g')
    %     % hold on 
    %     % plot(sigma70_strong_highConc1_empty_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_simData(:,species_oi_idx_list),2),...
    %     %     'LineWidth',3,'Color',[0.5,0.5,0.5])
    %     % plot(sigma70_strong_highConc1_empty_T7_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_T7_simData(:,species_oi_idx_list),2),...
    %     %     'LineWidth',2.25,'Color','r')
    %     % plot(sigma70_strong_highConc1_empty_sigma70_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_sigma70_simData(:,species_oi_idx_list),2),...
    %     %     'LineWidth',1.5,'Color','b')
    % 
    %     species_oi_plot_name = species_oi_plot_name_list{species_idx}; 
    %     title(species_oi_plot_name)
    %     xlabel('Time(hr)')
    %     ylabel('Concentration (nM)')
    %     set(gca,'FontSize',14)
    % 
    % 
    % end
    % legend('Baseline','Empty','Empty T7','Empty \sigma70')
    % 
    % % C. Resource competition 
    %     figure; 
    % 
    % species_oi_list_resource = {'RNAP','t7RNAP','Ribo','RNase','AGTP','CUTP','AA'};
    % for species_idx = 1:length(species_oi_list_resource)
    %     species_oi_name = species_oi_list_resource{species_idx};
    %     species_oi_idx_list = find(strcmp(species_name_list,species_oi_name));
    % 
    %     subplot(ceil(sqrt(length(species_oi_list_resource))),ceil(sqrt(length(species_oi_list_resource))),species_idx)
    %     % Plot low concentration 
    %     plot(sigma70_strong_lowConc2_no_empty_simTime ./ 3600,sum(sigma70_strong_lowConc2_no_empty_simData(:,species_oi_idx_list),2),...
    %         'LineWidth',3,'Color','g')
    %     hold on 
    %     plot(sigma70_strong_lowConc2_empty_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_simData(:,species_oi_idx_list),2),...
    %         'LineWidth',2.5,'Color',[0.5,0.5,0.5])
    %     plot(sigma70_strong_lowConc2_empty_T7_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_T7_simData(:,species_oi_idx_list),2),...
    %         'LineWidth',2,'Color','r')
    %     plot(sigma70_strong_lowConc2_empty_sigma70_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_sigma70_simData(:,species_oi_idx_list),2),...
    %         'LineWidth',1.5,'Color','b')
    % 
    %     % % Plot high concentration 
    %     % plot(sigma70_strong_highConc1_no_empty_simTime ./ 3600,sum(sigma70_strong_highConc1_no_empty_simData(:,species_oi_idx_list),2),'LineWidth',1.5,'Color','g')
    %     % hold on 
    %     % plot(sigma70_strong_highConc1_empty_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_simData(:,species_oi_idx_list),2),...
    %     %     'LineWidth',3,'Color',[0.5,0.5,0.5])
    %     % plot(sigma70_strong_highConc1_empty_T7_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_T7_simData(:,species_oi_idx_list),2),...
    %     %     'LineWidth',2.25,'Color','r')
    %     % plot(sigma70_strong_highConc1_empty_sigma70_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_sigma70_simData(:,species_oi_idx_list),2),...
    %     %     'LineWidth',1.5,'Color','b')
    % 
    % 
    %     title(species_oi_name)
    %     xlabel('Time(hr)')
    %     ylabel('Concentration (nM)')
    %     set(gca,'FontSize',14)
    % 
    % 
    % end
    % legend('Baseline','Empty','Empty T7','Empty sigma70')
    % 
    % % D. Breakdown of resource utilization 
    %     figure; 
    % 
    % species_oi_list_breakdown_resource = {'RNAP','RNase','Ribo'};
    % species_oi_list_breakdown_RNA = {'sfGFP','kanR','empty'}; 
    % num_plot = length(species_oi_list_breakdown_RNA) * length(species_oi_list_breakdown_resource);
    % species_oi_plot_name_list_bound_species = {'sfGFP DNA','kanR DNA','empty DNA',...
    %                                             'sfGFP mRNA','kanR mRNA','empty mRNA',...
    %                                             'sfGFP mRNA','kanR mRNA','empty mRNA'}; 
    % species_oi_plot_name_list_y_labels = {'RNAP bound','RNase bound','Ribosome bound'}; 
    % for species_idx_resource = 1:length(species_oi_list_breakdown_resource)
    %     species_oi_resource = species_oi_list_breakdown_resource{species_idx_resource};
    %     for species_idx_RNA = 1:length(species_oi_list_breakdown_RNA)
    %         species_oi_RNA = species_oi_list_breakdown_RNA{species_idx_RNA}; 
    %         plot_idx = (species_idx_resource - 1) * length(species_oi_list_breakdown_RNA) + species_idx_RNA; 
    %         subplot(ceil(sqrt(num_plot)),ceil(sqrt(num_plot)),plot_idx)
    % 
    %         if strcmp(species_oi_resource,'RNase')
    %             species_oi_idx_list = find(contains(species_name_list,species_oi_resource) &...
    %                 contains(species_name_list,species_oi_RNA));
    %         else
    %             species_oi_idx_list = find(contains(species_name_list,species_oi_resource) &...
    %                 contains(species_name_list,species_oi_RNA) & ~contains(species_name_list,'RNase'));
    %         end
    % 
    %         % Plot low concentration 
    %         plot(sigma70_strong_lowConc2_no_empty_simTime ./ 3600,sum(sigma70_strong_lowConc2_no_empty_simData(:,species_oi_idx_list),2),...
    %             'LineWidth',3,'Color','g')
    %         hold on 
    %         plot(sigma70_strong_lowConc2_empty_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_simData(:,species_oi_idx_list),2),...
    %             'LineWidth',2.5,'Color',[0.5,0.5,0.5])
    %         plot(sigma70_strong_lowConc2_empty_T7_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_T7_simData(:,species_oi_idx_list),2),...
    %             'LineWidth',2,'Color','r')
    %         plot(sigma70_strong_lowConc2_empty_sigma70_simTime ./ 3600,sum(sigma70_strong_lowConc2_empty_sigma70_simData(:,species_oi_idx_list),2),...
    %             'LineWidth',1.5,'Color','b')
    % 
    %         % % Plot high concentration 
    %         % plot(sigma70_strong_highConc1_no_empty_simTime ./ 3600,sum(sigma70_strong_highConc1_no_empty_simData(:,species_oi_idx_list),2),'LineWidth',1.5,'Color','g')
    %         % hold on 
    %         % plot(sigma70_strong_highConc1_empty_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_simData(:,species_oi_idx_list),2),...
    %         %     'LineWidth',3,'Color',[0.5,0.5,0.5])
    %         % plot(sigma70_strong_highConc1_empty_T7_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_T7_simData(:,species_oi_idx_list),2),...
    %         %     'LineWidth',2.25,'Color','r')
    %         % plot(sigma70_strong_highConc1_empty_sigma70_simTime ./ 3600,sum(sigma70_strong_highConc1_empty_sigma70_simData(:,species_oi_idx_list),2),...
    %         %     'LineWidth',1.5,'Color','b')
    % 
    %         species_oi_plot_name = species_oi_plot_name_list_bound_species{plot_idx}; 
    %         title(species_oi_plot_name)
    %         % if isequal(species_idx_RNA,1)
    %         %     ylabel_name = species_oi_plot_name_list_y_labels{species_idx_resource};
    %         % 
    %         % end
    %         xlabel('Time(hr)')
    %         ylabel('Concentration (nM)')
    %         set(gca,'FontSize',14)
    %     end
    % end
    % legend('Baseline','Empty','Empty T7','Empty sigma70')

%% Figure 4 
% 
%     % Relabel high resolution penalty labels 
%         % Label promotor + composition 
% relabel_list = cell(12,1); 
% for comp_idx = 1:length(composition_name_list)
%     composition_name = composition_name_list{comp_idx};
%     for prom_idx = 1:length(promotor_name_list)
%         promotor_name = promotor_name_list{prom_idx};
%         if contains(promotor_name,'sigma')
%             promotor_name_for_plot = strcat('\',strrep(promotor_name,'_',' '));
%         else
%             promotor_name_for_plot = strrep(promotor_name,'_',' ');
%         end
%         relabel_list{(comp_idx - 1) * length(promotor_name_list) + prom_idx} = ...
%             [promotor_name_for_plot ' ' composition_name ' ']; 
%     end
% end
%         % Add labels to high resolution label
% updated_high_res_penalty_labels = cell(32,1); 
% updated_high_res_penalty_labels(1:6) = high_res_penalty_term_labels(1:6);
% updated_high_res_penalty_labels(31:32) = high_res_penalty_term_labels(end-1:end); 
% positive_crosstalk_header = '(+) crosstalk ';
% positive_crosstalk_header_cell = cell(1,12);
% positive_crosstalk_header_cell(:) = {positive_crosstalk_header};
% negative_crosstalk_header = '(-) crosstalk ';
% negative_crosstalk_header_cell = cell(1,12);
% negative_crosstalk_header_cell(:) = {negative_crosstalk_header};
% updated_high_res_penalty_labels(7:18) = append(positive_crosstalk_header_cell',relabel_list); 
% updated_high_res_penalty_labels(19:30) = append(negative_crosstalk_header_cell',relabel_list); 
% 
% valid_high_res_penalty_terms = result_summary_file.all_high_res_penalty_terms(all(~isnan(result_summary_file.all_high_res_penalty_terms),2),:); 
% valid_all_sampled_params = init_param_file.satisfied_sampled_params(all(~isnan(result_summary_file.all_high_res_penalty_terms),2),:);
% rank_valid_all_penalty_terms = tiedrank(valid_high_res_penalty_terms);
% rank_param_mat = tiedrank(valid_all_sampled_params); 
% augmented_PRCC = partialcorri(rank_valid_all_penalty_terms,rank_param_mat); 
% 
% figure;                                                                                                                                                                                                                                                            
% h3 = heatmap(strrep(parameters_name_list,'_',' '),strrep(updated_high_res_penalty_labels,'_',' '),augmented_PRCC,'Colormap',jet); 
% xlabel('Parameter Names')
% ylabel('Qualitative Terms') 
% title('Partial Rank Correlation Coefficient')
% set(gca,'FontSize',14)
%% Figure 5 

    % Result file 
TX_gen_on_RNAP_resultFileName ='result_files/20231101_param_est_run0929_RNAP_toxin_mech.mat';
TX_gen_on_Ribo_resultFileName = 'result_files/20230125_param_est_run1111_TXgen_Ribo.mat';
TL_gen_on_RNAP_resultFileName = 'result_files/20230125_param_est_run2230_TLgen_RNAP.mat';
TL_gen_on_Ribo_resultFileName = 'result_files/20230125_param_est_run1326_TLgen_Ribo.mat';

TX_gen_on_RNAP_resultFile = load(TX_gen_on_RNAP_resultFileName);
TX_gen_on_Ribo_resultFile = load(TX_gen_on_Ribo_resultFileName); 
TL_gen_on_RNAP_resultFile = load(TL_gen_on_RNAP_resultFileName);
TL_gen_on_Ribo_resultFile = load(TL_gen_on_Ribo_resultFileName);

    % TX generated toxin poisons RNAP
TX_gen_on_RNAP_metric_summary = TX_gen_on_RNAP_resultFile.metric_summary; 
[~,max_idx] = max(TX_gen_on_RNAP_metric_summary); 
TX_gen_on_RNAP_fitResult_oi = TX_gen_on_RNAP_resultFile.all_fitResults{max_idx,1};
fitted_data_TX_gen_on_RNAP = fitted(TX_gen_on_RNAP_fitResult_oi); 
[Time_TX_gen_on_RNAP,Data_TX_gen_on_RNAP] = process_crosstalk_data_from_source(fitted_data_TX_gen_on_RNAP,'simulated_data');

    % TX generated toxin poisons Ribosome 
TX_gen_on_Ribo_metric_summary = TX_gen_on_Ribo_resultFile.metric_summary; 
[~,max_idx] = max(TX_gen_on_Ribo_metric_summary); 
TX_gen_on_Ribo_fitResult_oi = TX_gen_on_Ribo_resultFile.all_fitResults{max_idx,1};
fitted_data_TX_gen_on_Ribo = fitted(TX_gen_on_Ribo_fitResult_oi); 
[Time_TX_gen_on_Ribo,Data_TX_gen_on_Ribo] = process_crosstalk_data_from_source(fitted_data_TX_gen_on_Ribo,'simulated_data');
        %%%% Temp diagnosis %%%%
        figure; 
        for idx = 1:length(TX_gen_on_Ribo_resultFile.all_fitResults)
            fitResult_oi = TX_gen_on_Ribo_resultFile.all_fitResults{idx};
            fitted_data = fitted(fitResult_oi); 
            [Time,Data] = process_crosstalk_data_from_source(fitted_data,'simulated_data');
            subplot(7,7,idx)
            for conc_idx = 1:length(conc_vec)
                selected_time = Time{conc_idx};
                selected_data = Data{conc_idx}; 
                plot(selected_time,selected_data,'LineWidth',1.5)
                hold on 
            end

        end






        TX_gen_on_Ribo_problemObject = TX_gen_on_Ribo_resultFile.problemObject;
        Model_TX_gen_on_Ribo = TX_gen_on_Ribo_problemObject.Model;
        species_name_list = {Model_TX_gen_on_Ribo.Species.Name};
        species_oi_name = 'protein sfGFP*';
        species_oi_idx = find(strcmp(species_name_list,species_oi_name)); 
        additional_track_species = [species_name_list(1:species_oi_idx -1) species_name_list(species_oi_idx + 1:end)];
        simFunction_TX_gen_on_Ribo = create_simFun_from_problemObject(TX_gen_on_Ribo_problemObject,additional_track_species); 
        dosing_information = create_dosing_info_from_problemObject(TX_gen_on_Ribo_problemObject); 
        estimated_params = [TX_gen_on_Ribo_fitResult_oi.ParameterEstimates.Estimate]';
        [Time,Data] = simFunction_TX_gen_on_Ribo(estimated_params,tEnd,dosing_information,tStart:tEnd); 

        simFunction_species_name_list = {simFunction_TX_gen_on_Ribo.Observables.Name}; 
        simFunction_species_name_list = simFunction_species_name_list{1};
        conc_idx_list = [1,3,7];
        figure; 
        for conc_idx = conc_idx_list
        selected_time = Time{conc_idx};
        selected_data = Data{conc_idx}; 
        species_oi_name_list = {'protein sfGFP*','protein sfGFP','RNA utrGFP--sfGFP','toxin','Ribo','RNase'}; 
        species_oi_name_for_plot_list = {'mature sfGFP','unmatured sfGFP','sfGFP mRNA','Toxin','Ribosome','RNase'};
        % figure; 
        % for species_oi_idx = 1:length(simFunction_species_name_list)
        %     species_oi_name = simFunction_species_name_list{species_oi_idx};
        %     % species_idx = strcmp(simFunction_species_name_list,species_oi_name);
        %     subplot(10,10,species_oi_idx)
        %     plot(selected_time ./ 3600,selected_data(:,species_oi_idx),'LineWidth',2)
        %     hold on
        %     xlabel('Time (hr)')
        %     ylabel('Concentration (nM)')
        %     title(strrep(simFunction_species_name_list{species_oi_idx},'_',' '))
        %     set(gca,'FontSize',14)
        % end

        for species_oi_idx = 1:length(species_oi_name_list)
            species_oi_name = species_oi_name_list{species_oi_idx};
            species_idx = strcmp(simFunction_species_name_list,species_oi_name);
            subplot(2,3,species_oi_idx)
            plot(selected_time ./ 3600,selected_data(:,species_idx),'LineWidth',2)
            hold on
            xlabel('Time (hr)')
            ylabel('Concentration (nM)')
            title(species_oi_name_for_plot_list{species_oi_idx})
            set(gca,'FontSize',14)
        end
        end
        %%%%% end of temp diagnosis %%%%%
        

    % TL generated toxin poisons RNAP
TL_gen_on_RNAP_metric_summary = TL_gen_on_RNAP_resultFile.metric_summary; 
[~,max_idx] = max(TL_gen_on_RNAP_metric_summary); 
TL_gen_on_RNAP_fitResult_oi = TL_gen_on_RNAP_resultFile.all_fitResults{max_idx,1};
fitted_data_TL_gen_on_RNAP = fitted(TL_gen_on_RNAP_fitResult_oi); 
[Time_TL_gen_on_RNAP,Data_TL_gen_on_RNAP] = process_crosstalk_data_from_source(fitted_data_TL_gen_on_RNAP,'simulated_data');

    % TL generated toxin poisons Ribosome 
TL_gen_on_Ribo_metric_summary = TL_gen_on_Ribo_resultFile.metric_summary; 
[~,max_idx] = max(TL_gen_on_Ribo_metric_summary); 
TL_gen_on_Ribo_fitResult_oi = TL_gen_on_Ribo_resultFile.all_fitResults{max_idx,1};
fitted_data_TL_gen_on_Ribo = fitted(TL_gen_on_Ribo_fitResult_oi); 
[Time_TL_gen_on_Ribo,Data_TL_gen_on_Ribo] = process_crosstalk_data_from_source(fitted_data_TL_gen_on_Ribo,'simulated_data');

% figure; 
% subplot(2,2,1)
% for conc_idx = 1:length(Data_TX_gen_on_RNAP)
%     plot(Time_TX_gen_on_RNAP{conc_idx}./3600,Data_TX_gen_on_RNAP{conc_idx},'LineWidth',1.5,'Color',sampledColors(conc_idx,:))
%     hold on 
% end
% xlabel('Time (hr)','FontSize',14)
% ylabel('Concentration (nM)','FontSize',14)
% set(gca,'FontSize',14)
% 
% subplot(2,2,2)
% for conc_idx = 1:length(Data_TL_gen_on_Ribo)
%     plot(Time_TL_gen_on_Ribo{conc_idx}./3600,Data_TL_gen_on_Ribo{conc_idx},'LineWidth',1.5,'Color',sampledColors(conc_idx,:))
%     hold on 
% end
% xlabel('Time (hr)','FontSize',14)
% ylabel('Concentration (nM)','FontSize',14)
% set(gca,'FontSize',14)
% 
% subplot(2,2,3)
% for conc_idx = 1:length(Data_TL_gen_on_RNAP)
%     plot(Time_TL_gen_on_RNAP{conc_idx}./3600,Data_TL_gen_on_RNAP{conc_idx},'LineWidth',1.5,'Color',sampledColors(conc_idx,:))
%     hold on 
% end
% xlabel('Time (hr)','FontSize',14)
% ylabel('Concentration (nM)','FontSize',14)
% set(gca,'FontSize',14)
% 
% 
% subplot(2,2,4)
% for conc_idx = 1:length(Data_TX_gen_on_Ribo)
%     plot(Time_TX_gen_on_Ribo{conc_idx}./3600,Data_TX_gen_on_Ribo{conc_idx},'LineWidth',1.5,'Color',sampledColors(conc_idx,:))
%     hold on 
% end
% xlabel('Time (hr)','FontSize',14)
% ylabel('Concentration (nM)','FontSize',14)
% set(gca,'FontSize',14)


%% Figure 6 

% high_res_penalty_corr_Spearman = corr(valid_high_res_penalty_terms,'Type','Spearman'); 
% high_res_penalty_corr = corr(valid_high_res_penalty_terms); 
% 
% figure; 
% h3 = heatmap(strrep(updated_high_res_penalty_labels,'_',' '),strrep(updated_high_res_penalty_labels,'_',' '),...
%     high_res_penalty_corr_Spearman,'Colormap',jet); 
% xlabel('Penalty Terms')
% ylabel('Penalty Terms') 
% title('Spearman Correlation Coefficient')

%% Figure S3
%     % Plot baseline and crosstalk ratio for experimental data 
% experimental_data = problemObject.Data;
% [simulated_time,simulated_data] = process_crosstalk_data_from_source(experimental_data,'grouped_data');
% plot_simulated(simulated_time,simulated_data,'baseline_short')
% plot_simulated(simulated_time,simulated_data,'crosstalk_ratio')


%% Figure S4, S10, S13
% 
%     % Figure S4 - original model 
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
% plot_qual_obj_dist(all_high_res_penalty_term,penalty_cutoff)

%     % Figure S10 - satisfy at least 8 
% satisfy_flag_matrix = streamline_penalty_criteria(result_summary_file.all_high_res_penalty_term,penalty_cutoff);
% num_satisfied_criteria_per_sample = sum(satisfy_flag_matrix,2); 
% num_satisfied_8_sample_idx_list = find(num_satisfied_criteria_per_sample >= 8); 
% high_res_penalty_terms_oi = result_summary_file.all_high_res_penalty_term(num_satisfied_8_sample_idx_list,:); 
% plot_qual_obj_dist(high_res_penalty_terms_oi,penalty_cutoff)

%     % Figure S13 - surrogate term as obj function 
% result_file_qual_obj = load('test_save_files/20240919_lhd_sampled_params_optimizatio_qual_obj_result_summary.mat'); 
% satisfy_flag_matrix_qual_obj = streamline_penalty_criteria(result_file_qual_obj.all_high_res_penalty_qual_obj,penalty_cutoff);
% 
% num_satisfied_criteria_per_sample_streamlined_qual_obj = sum(satisfy_flag_matrix_qual_obj,2); 
% filter_idx_list_qual_obj = find(num_satisfied_criteria_per_sample_streamlined_qual_obj >= 8); 
% filtered_high_res_penalty_qual_obj = result_file_qual_obj.all_high_res_penalty_qual_obj(filter_idx_list_qual_obj,:);
% 
% plot_qual_obj_dist(filtered_high_res_penalty_qual_obj,penalty_cutoff)

%% Figure S6

% RNA_deg_fitting_result = load('test_save_files/20240916_alternative_kinetics_fitting_original_data_truobleshooted_mechanism.mat'); 
% RNA_deg_fitting_result_part2 = load('test_save_files/20240916_alternative_kinetics_fitting_original_data_truobleshooted_mechanism_part2.mat'); 
% RNA_deg_fitting_result_part3 = load('test_save_files/20240916_alternative_kinetics_fitting_original_data_truobleshooted_mechanism_part3.mat'); 
% data_for_plot = [min(RNA_deg_fitting_result_part2.all_fval_firstOrder),...
%     min(RNA_deg_fitting_result_part2.all_fval_MichaelisMenten),...
%     min(RNA_deg_fitting_result_part2.all_fval_MassAction),...
%     min(RNA_deg_fitting_result.all_fval_TestPrototype),...
%     min(RNA_deg_fitting_result_part3.all_fval_TestBindingSites)];
% xlabels = {'1st order','Michaelis-Menten','Mass Action','Further Degradation',...
%     'Binding Site'}; 
% xlabels_for_plot = categorical(xlabels);
% xlabels_for_plot = reordercats(xlabels_for_plot,xlabels); 
% figure; 
% bar(xlabels_for_plot,data_for_plot)
% ylabel('minimum SSE');
% set(gca,'FontSize',14)
%% Figure S7 & S8
% result_file_name_wo_weights = 'result_files/20240829_param_est_run1207_wo_weight.mat';
% result_file = load(result_file_name_wo_weights);
% 
% all_SSE = result_file.all_SSE;
% all_high_res_penalty_terms = result_file.all_high_res_penalty_terms; 
% all_init_SSE = result_file.all_init_SSE;
% all_init_penalty_term = result_file.all_init_penalty_term; 
% 
% satisfy_flag_matrix_w_optimization = streamline_penalty_criteria(all_high_res_penalty_terms,penalty_cutoff);
% satisfy_flag_matrix_wo_optimization = streamline_penalty_criteria(all_init_penalty_term,penalty_cutoff);
% 
% figure; 
% scatter(log10(all_init_SSE),sum(satisfy_flag_matrix_wo_optimization,2),'filled')
% hold on 
% scatter(log10(all_SSE),sum(satisfy_flag_matrix_w_optimization,2),'filled')
% xlabel('log10(SSE)')
% ylabel('# satisfied')
% legend('Sampled','Optimized')
% set(gca,'FontSize',14)
% 
% result_file_name_w_weights = 'param_est_run_save/20240829_param_est_run1207_w_weight.mat';
% result_file_w_weights = load(rsetesult_file_name_w_weights); 
% all_high_res_penalty_terms_w_weights = result_file_w_weights.all_high_res_penalty_terms; 
% 
% satisfy_flag_matrix_wo_weights = streamline_penalty_criteria(all_high_res_penalty_terms,penalty_cutoff);
% satisfy_flag_matrix_w_weights = streamline_penalty_criteria(all_high_res_penalty_terms_w_weights,penalty_cutoff);
% 
% % Compile mastrix for bar plot 
% num_satisfied_for_plot = nan(size(satisfy_flag_matrix_wo_weights,2),2);
% for criteria_idx = 1:size(satisfy_flag_matrix_wo_weights,2)
%     num_satisfied_for_plot(criteria_idx,1) = sum(satisfy_flag_matrix_wo_weights(:,criteria_idx)); 
%     num_satisfied_for_plot(criteria_idx,2) = sum(satisfy_flag_matrix_w_weights(:,criteria_idx)); 
% 
% end
% 
% figure; 
% xlabels = categorical(criteria_name_list); 
% xlabels = reordercats(xlabels,criteria_name_list); 
% bar(xlabels,num_satisfied_for_plot)
% legend('Control','Normalized data')
% set(gca,'FontSize',14)
% ylabel('# satisfied')
%% Figure S9
    % Funnel chart 

    % Calculate number of criteria each row meets
    % satisfy_flag_matrix = streamline_penalty_criteria(result_summary_file.all_high_res_penalty_terms,penalty_cutoff);
    % 
    % num_satisfied_criteria_per_sample = sum(satisfy_flag_matrix,2); 
    % num_satisfied_sample_all = nan(size(satisfy_flag_matrix,2),1);
    % for num_satisfied = 1:size(satisfy_flag_matrix,2)
    %     num_satisfied_sample = sum(num_satisfied_criteria_per_sample >= num_satisfied); 
    %     num_satisfied_sample_all(num_satisfied) = num_satisfied_sample; 
    % end

    % funnelchart(num_satisfied_sample_all)
    % ylabel('Number of criteria satisfied')
    % xlabel('Number of samples')

%% Figure S11
%     % Plot out TX crosstalk for selected model 
% TX_group_description = {'T7_strong_no_empty','T7_strong_empty','T7_strong_empty_T7','T7_strong_empty_sigma70',...
%     'T7_weak_no_empty','T7_weak_empty','T7_weak_empty_T7','T7_weak_empty_sigma70',...
%     'sigma70_strong_no_empty','sigma70_strong_empty','sigma70_strong_empty_T7','sigma70_strong_empty_sigma70'};
% param_info_path = 'parameters_v2_expanded.xlsx';
% problemObject_TX = setProblemObject_v3(TX_group_description,[],param_info_path);
% dosing_information_TX = create_dosing_info_from_problemObject(problemObject_TX); 
% 
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

% 
%     % Load fitted TX-level crosstalk results 
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
% end
% 
% [~,min_idx_2] = min(metric_summary_updated); 
% opt_fitResult_updated = TX_fitted_result_file_updated.all_fitResults{min_idx_2}; 
% simData_updated = fitted(opt_fitResult_updated);
% [simulated_time_updated,simulated_data_updated] = process_crosstalk_data_from_source(simData_updated,'simulated_data');
% 
% plot_time_course_crosstalk(simulated_time_updated,simulated_data_updated)





