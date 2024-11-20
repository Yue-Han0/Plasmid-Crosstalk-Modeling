clear
clc

% This script contains code to fit the mRNA degradtion data using
% alternative mechanisms (1) 1st-order (2) MM (3) mass-action (4) Mixed
% processive and non-processive degradation 

currentpath = pwd; 
addpath(genpath(currentpath))
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))

%% Define Simbiology Models assuming different kinetics 

all_kinetics = {'firstOrder','MichaelisMenten','MassAction','MixedDeg','TestPrototype','TestBindingSites'}; 
for kinetics_idx = 1:length(all_kinetics)

    kinetics = all_kinetics{kinetics_idx}; 
    eval(sprintf('Model_%s = build_RNAdeg_SimBiology_Model(kinetics);',kinetics)); 

end

%% Load experimental data 
gain_70_timeVec = readtable('processed_data/degradation_curve_gain_70.xlsx','Sheet','timeVec');
gain_70_RNA_time_course = readtable('processed_data/degradation_curve_gain_70.xlsx','Sheet','degradation');
gain_75_timeVec = readtable('processed_data/degradation_curve_gain_75.xlsx','Sheet','timeVec');
gain_75_RNA_time_course = readtable('processed_data/degradation_curve_gain_75.xlsx','Sheet','degradation');

        % Get initial concentration 
RNA_init_conc_gain70 = gain_70_RNA_time_course{1,2:9}; 
RNA_init_conc_gain75 = gain_75_RNA_time_course{1,2:8}; 
        % Get timeVec and time-course 
RNA_timeVec_gain70 = gain_70_timeVec{:,2}; 
RNA_timeVec_gain75 = gain_75_timeVec{:,2}; 
RNA_fluo_time_course_gain70 = gain_70_RNA_time_course{2:end,11:18}; 
RNA_fluo_time_course_gain75 = gain_75_RNA_time_course{2:end,10:16}; 

   % From imputed processed data 
% imputed_processed_RNA_deg_file = load('imputed_RNAdeg_data.mat'); 
%     % Use imputed data 
% RNA_fluo_time_course_gain70 = imputed_processed_RNA_deg_file.imputed_data_gain70; 
% RNA_init_conc_gain70 = imputed_processed_RNA_deg_file.imputed_init_cond_gain70; 
% RNA_timeVec_gain70 = imputed_processed_RNA_deg_file.gain_70_deg_data_timeVec{:,2}; 
    % Use imputed simulated data (all simulation) 
% RNA_fluo_time_course_gain70 = imputed_processed_RNA_deg_file.sim_imputed_data_gain70; 
% RNA_init_conc_gain70  = imputed_processed_RNA_deg_file.sim_imputed_init_cond_gain70;
% RNA_timeVec_gain70 = imputed_processed_RNA_deg_file.full_sim_timeVec_gain70'; 

        % Convert fluorescence to mRNA concentration based on calibration
        % curve
% Note: these two should technically give the same mRNA concentrations, but
% due to experimental noise and fitting erros there's a small difference
% (~10%) between them. Accounting for uncertainty propogation here would be
% insane, so let's use gain 70 to fit degradation. 
RNA_fluo_time_course_gain70 = RNA_fluo_time_course_gain70(:,~any(isnan(RNA_fluo_time_course_gain70),1));
RNA_time_course_gain70 = fluo_mRNA_conversion(RNA_fluo_time_course_gain70,[],70); 

    % % Add a time delay to timeVec and data
RNA_timeVec_gain70 = RNA_timeVec_gain70 + 4 * 60; % 4-min time delay 
RNA_timeVec_gain70 = [0;RNA_timeVec_gain70]; % Add a true time zero 
RNA_time_course_gain70 = [RNA_init_conc_gain70;RNA_time_course_gain70]; % Add init conc for the true time zero 

    % Process data into a Simbiology readable format 
Time = repmat(RNA_timeVec_gain70,length(RNA_init_conc_gain70),1);
mRNA_concentration = RNA_time_course_gain70(:); 
    % Treat negative values as 0 
mRNA_concentration(mRNA_concentration < 0) = 0; 

Group = nan(length(mRNA_concentration),1); 
mRNA_concentration_add = zeros(length(mRNA_concentration),1); 
for group_idx = 1:length(RNA_init_conc_gain70)
    Group((group_idx - 1) * length(RNA_timeVec_gain70) + 1:group_idx * length(RNA_timeVec_gain70)) = group_idx; 
    mRNA_concentration_add((group_idx - 1) * length(RNA_timeVec_gain70) + 1) = RNA_init_conc_gain70(group_idx); 
end
RNA_degradation_data_table = array2table([Group,Time,mRNA_concentration,mRNA_concentration_add], 'VariableNames', {'Group','Time', 'mRNA_concentration','mRNA_concentration_add'});

%% Create a fitProblem object 
for kinetics_idx = 1:length(all_kinetics)

    kinetics = all_kinetics{kinetics_idx}; 
    eval(sprintf('problemObject_%s = build_RNAdeg_fitProblem_struct(kinetics,Model_%s,RNA_degradation_data_table);',kinetics,kinetics)); 

end

%% Fit using simFunction 
% for kinetics_idx = length(all_kinetics):length(all_kinetics)
% 
%     kinetics = all_kinetics{kinetics_idx}; 
%     eval(sprintf('[all_fval_%s,all_estimated_params_%s] = fit_RNAdeg_exp_data(Model_%s,problemObject_%s);',kinetics,kinetics,kinetics,kinetics)); 
% 
% end


% save('test_save_files/20240916_alternative_kinetics_fitting_original_data_truobleshooted_mechanism_part3.mat'); 
% save('test_save_files/20240912_alterantive_kinetics_fitting_original_data_expanded_bound.mat'); 
%% Analysis 
% load('test_save_files/20240912_alterantive_kinetics_fitting_imputed_data_expanded_bound.mat'); 
load('test_save_files/20240916_alternative_kinetics_fitting_original_data_truobleshooted_mechanism.mat');
load('test_save_files/20240916_alternative_kinetics_fitting_original_data_truobleshooted_mechanism_part2.mat');
load('test_save_files/20240916_alternative_kinetics_fitting_original_data_truobleshooted_mechanism_part3.mat'); 
    %%  Plot out SSE distirbution 

figure;
for kinetics_idx = 1:length(all_kinetics)
    subplot(3,2,kinetics_idx)
    kinetics = all_kinetics{kinetics_idx}; 
    eval(sprintf('histogram(log10(all_fval_%s))',kinetics))
    xlabel('log_{10}SSE')
    ylabel('# Occurence')
    title(kinetics)
    xlim([7.5,11])
end

    %% Compare the minimal SSE 
min_SSEs = [min(all_fval_firstOrder),min(all_fval_MichaelisMenten),min(all_fval_MassAction),min(all_fval_MixedDeg),min(all_fval_TestPrototype),min(all_fval_TestBindingSites)];
figure;
xlabels = {'First Order','Michaelis-Menten','Mass Action','New Mechanism v1','New Mechanism v2','Binding Site Mechanism'}; 
xlabels_for_plot = categorical(xlabels);
xlabels_for_plot = reordercats(xlabels_for_plot,xlabels); 
bar(xlabels_for_plot,min_SSEs)


    %% Experimental vs. simulated & estimated parameter distribution 
% dosing_information = create_dosing_info_from_problemObject(problemObject_firstOrder); 
% tStart = 0; 
% tEnd = 2940; 
% kinetics_logSSE_cutoff_list = [1e+08,1e+08,1e+08,1e+08,1e+08,1e+08]; 
% for kinetic_idx = 1:length(all_kinetics)
% 
%     kinetics = all_kinetics{kinetic_idx}; 
% 
%         % Plot out simulated vs. experimental
%     eval(sprintf('[~,sort_idx] = sort(all_fval_%s,''ascend'');',kinetics))
%     eval(sprintf('opt_estimated_params = all_estimated_params_%s(sort_idx(1),:);',kinetics))
%     eval(sprintf('simFunction = create_simFun_from_problemObject(problemObject_%s);',kinetics));
%     [simTime,simData] = simFunction(opt_estimated_params,tEnd,dosing_information,tStart:tEnd); 
%     figure; 
%     for conc_idx = 1:length(simData)
%         subplot(ceil(sqrt(length(simData))),ceil(sqrt(length(simData))),conc_idx)
%         plot(simTime{conc_idx},simData{conc_idx}(:,1),'LineWidth',1.5,'Color','r')
%         hold on 
%         plot(RNA_degradation_data_table.Time(RNA_degradation_data_table.Group == conc_idx),...
%             RNA_degradation_data_table.mRNA_concentration(RNA_degradation_data_table.Group == conc_idx),'LineWidth',1.5,'Color','k')
%     end
%     sgtitle(kinetics)
%     % saveas(gcf,sprintf('plots/20240807_RNAdegFit_SimvExp_%s.png',kinetics),'png')
% 
%         % Check parameter distribution to assess parameter sloppiness 
%     eval(sprintf('all_fval = all_fval_%s;',kinetics));
%     eval(sprintf('all_estimated_params = all_estimated_params_%s;',kinetics)); 
%     eval(sprintf('RNA_deg_probObject = problemObject_%s;',kinetics)); 
%     figure; 
%     kinetics_logSSE_cutoff = kinetics_logSSE_cutoff_list(kinetic_idx); 
%     lower_SSE_idx_list = (all_fval < kinetics_logSSE_cutoff); 
%     lower_SSE_params = all_estimated_params(lower_SSE_idx_list,:); 
%     higher_SSE_params = all_estimated_params(~lower_SSE_idx_list,:); 
%     param_names = {RNA_deg_probObject.Estimated.Name};
%     for param_idx = 1:length(RNA_deg_probObject.Estimated)
% 
%         subplot(ceil(sqrt(length(RNA_deg_probObject.Estimated))),ceil(sqrt(length(RNA_deg_probObject.Estimated))),param_idx)
%         histogram(log(lower_SSE_params(:,param_idx)))
%         hold on
%         % histogram(log(higher_SSE_params(:,param_idx)))
%         title(strrep(param_names{param_idx},'_',' '))
% 
%     end
%     sgtitle(kinetics)
%     % legend('Lower SSE','Higher SSE')
%     % saveas(gcf,sprintf('plots/20240807_RNAdegFit_ParamDist_%s.png',kinetics),'png')
% 
% end

%%  Detailed time-course of species concentrations and reaction rates 

for kinetics_idx = 6:length(all_kinetics)

    kinetics = all_kinetics{kinetics_idx}; 

        % Take the optimal fit
    eval(sprintf('[~,opt_idx_%s] = min(all_fval_%s);',kinetics,kinetics));
    eval(sprintf('opt_estimated_params_oi = all_estimated_params_%s(opt_idx_%s,:);',kinetics,kinetics)); 

        % Create simFunction & Run simulation with optimal parameters 
    eval(sprintf('model_species = {Model_%s.Species.Name}; ',kinetics));
    species_oi_name = 'RNA utrbroc--no_protein'; 
    species_oi_idx = find(strcmp(species_oi_name,model_species)); 
    additional_track_species = [model_species(1:species_oi_idx - 1),model_species(species_oi_idx + 1:end)];
    eval(sprintf('problemObject_oi = problemObject_%s;',kinetics)); 
    simFunction = create_simFun_from_problemObject(problemObject_oi,additional_track_species); 
    dosing_information = create_dosing_info_from_problemObject(problemObject_oi);

    tEnd = 3000; 
    tStart = 0; 

        % Check optimal estimated params 
    [simTime,simData] = simFunction(opt_estimated_params_oi,tEnd,dosing_information,tStart:tEnd);
        % Take the 1st mRNA concentration as case study 

    % % In this optimal case, mRNA was not degraded since TXTL_RNAdeg_kc_proc is
    % % not larger than the other kinetic parameters 
    %     % Let's look at cases where TXTL_RNAdeg_kc_proc is larger than other
    %     % kinetic parameters
    % select_idx_list = find(all_estimated_params_MixedDeg(:,4) > all_estimated_params_MixedDeg(:,1) &...
    %     all_estimated_params_MixedDeg(:,4) > all_estimated_params_MixedDeg(:,2));
    % selected_estimated_params = all_estimated_params_MixedDeg(select_idx_list,:); 
    % [simTime,simData] = simFunction(selected_estimated_params(1,:),tEnd,dosing_information,tStart:tEnd);
        % Take the 1st mRNA concentration as case study 

        % Check out customized parameters 
    % customized_params = [0.01,0.01,100,0.01];
    % [simTime,simData] = simFunction(customized_params,tEnd,dosing_information,tStart:tEnd);

    simTime_lowestConc = simTime{1};
    simData_lowestConc = simData{1}; 
    simTime_highestConc = simTime{end};
    simData_highestConc = simData{end}; 

    master_title_prefix = sprintf('%s kinetics',kinetics);

        % (1) plot out all species time-course for lowest conc 
    figure; 
    all_tracked_species = [{species_oi_name},additional_track_species]; 
    for species_idx = 1:length(all_tracked_species)
        tracked_species_name = all_tracked_species{species_idx};
        species_time_course = simData_lowestConc(:,species_idx);

        subplot(ceil(sqrt(length(all_tracked_species))),ceil(sqrt(length(all_tracked_species))),species_idx)
        plot(simTime_lowestConc(2:end),species_time_course(2:end),'LineWidth',1.5)
        title(strrep(tracked_species_name,'_',' '))

    end
    sgtitle(sprintf('%s Full Time-Course Lowest mRNA concentration',master_title_prefix)); 

            % (2) plot out all species time-course for lowest conc -
            % truncated time-course
    figure; 
    all_tracked_species = [{species_oi_name},additional_track_species]; 
    for species_idx = 1:length(all_tracked_species)
        tracked_species_name = all_tracked_species{species_idx};
        species_time_course = simData_lowestConc(:,species_idx);

        subplot(ceil(sqrt(length(all_tracked_species))),ceil(sqrt(length(all_tracked_species))),species_idx)
        truncate_idx = find(simTime_lowestConc > 1); 
        plot(simTime_lowestConc(truncate_idx(1):end),species_time_course(truncate_idx(1):end),'LineWidth',1.5)
        title(strrep(tracked_species_name,'_',' '))

    end
    sgtitle(sprintf('%s Truncated Time-Course Lowest mRNA concentration',master_title_prefix)); 

            % (3) plot out all species time-course for highest conc 
    figure; 
    all_tracked_species = [{species_oi_name},additional_track_species]; 
    for species_idx = 1:length(all_tracked_species)
        tracked_species_name = all_tracked_species{species_idx};
        species_time_course = simData_highestConc(:,species_idx);

        subplot(ceil(sqrt(length(all_tracked_species))),ceil(sqrt(length(all_tracked_species))),species_idx)
        plot(simTime_highestConc(2:end),species_time_course(2:end),'LineWidth',1.5)
        title(strrep(tracked_species_name,'_',' '))

    end
    sgtitle(sprintf('%s Full Time-Course Highest mRNA concentration',master_title_prefix)); 

            % (4) plot out all species time-course for highest conc -
            % truncated time-course
    figure; 
    all_tracked_species = [{species_oi_name},additional_track_species]; 
    for species_idx = 1:length(all_tracked_species)
        tracked_species_name = all_tracked_species{species_idx};
        species_time_course = simData_highestConc(:,species_idx);

        subplot(ceil(sqrt(length(all_tracked_species))),ceil(sqrt(length(all_tracked_species))),species_idx)
        truncate_idx = find(simTime_highestConc > 1); 
        plot(simTime_highestConc(truncate_idx(1):end),species_time_course(truncate_idx(1):end),'LineWidth',1.5)
        title(strrep(tracked_species_name,'_',' '))

    end
    sgtitle(sprintf('%s Truncated Time-Course Highest mRNA concentration',master_title_prefix)); 

        % Calculate reaction rate time-course as well 
    estimated_parameter_name_list = {problemObject_oi.Estimated.Name};
    eval(sprintf('[rxn_name_list,rxn_time_course_lowConc] = get_reaction_time_course(Model_%s,estimated_parameter_name_list,opt_estimated_params_oi,simTime_lowestConc,simData_lowestConc,all_tracked_species);',kinetics)); 
    eval(sprintf('[~,rxn_time_course_highConc] = get_reaction_time_course(Model_%s,estimated_parameter_name_list,opt_estimated_params_oi,simTime_highestConc,simData_highestConc,all_tracked_species);',kinetics)); 

    % (1) Plot out all reaction full time-course for lowest conc
    figure; 
    for rxn_idx = 1:length(rxn_name_list)
        rxn_name = rxn_name_list{rxn_idx};
        rxn_time_course_single = rxn_time_course_lowConc(:,rxn_idx); 

        subplot(ceil(sqrt(length(rxn_name_list))),ceil(sqrt(length(rxn_name_list))),rxn_idx)
        plot(simTime_lowestConc(2:end),rxn_time_course_single(2:end),'LineWidth',1.5)
        title(strrep(rxn_name,'_',' '))

    end
    sgtitle(sprintf('%s Reaction Full time-course Lowest mRNA concentration',master_title_prefix))

        % (2) Plot out all reaction truncated time-course for lowest conc
    figure; 
    for rxn_idx = 1:length(rxn_name_list)
        rxn_name = rxn_name_list{rxn_idx};
        rxn_time_course_single = rxn_time_course_lowConc(:,rxn_idx); 

        subplot(ceil(sqrt(length(rxn_name_list))),ceil(sqrt(length(rxn_name_list))),rxn_idx)
        truncate_idx = find(simTime_lowestConc > 1); 
        plot(simTime_lowestConc(truncate_idx(1):end),rxn_time_course_single(truncate_idx(1):end),'LineWidth',1.5)
        title(strrep(rxn_name,'_',' '))

    end
    sgtitle(sprintf('%s Reaction Truncated time-course Lowest mRNA concentration',master_title_prefix))

        % (3) Plot out all reaction full time-course for highest conc
    figure; 
    for rxn_idx = 1:length(rxn_name_list)
        rxn_name = rxn_name_list{rxn_idx};
        rxn_time_course_single = rxn_time_course_highConc(:,rxn_idx); 

        subplot(ceil(sqrt(length(rxn_name_list))),ceil(sqrt(length(rxn_name_list))),rxn_idx)
        plot(simTime_highestConc(2:end),rxn_time_course_single(2:end),'LineWidth',1.5)
        title(strrep(rxn_name,'_',' '))

    end
    sgtitle(sprintf('%s Reaction Full time-course Highest mRNA concentration',master_title_prefix))

        % (4) Plot out all reaction truncated time-course for highest conc
    figure; 
    for rxn_idx = 1:length(rxn_name_list)
        rxn_name = rxn_name_list{rxn_idx};
        rxn_time_course_single = rxn_time_course_highConc(:,rxn_idx); 

        subplot(ceil(sqrt(length(rxn_name_list))),ceil(sqrt(length(rxn_name_list))),rxn_idx)
        truncate_idx = find(simTime_highestConc > 1); 
        plot(simTime_highestConc(truncate_idx(1):end),rxn_time_course_single(truncate_idx(1):end),'LineWidth',1.5)
        title(strrep(rxn_name,'_',' '))

    end
    sgtitle(sprintf('%s Reaction Truncated time-course Highest mRNA concentration',master_title_prefix))


end
