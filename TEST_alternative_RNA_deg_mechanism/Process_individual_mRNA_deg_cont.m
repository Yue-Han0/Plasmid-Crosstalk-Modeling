clear
clc

% addpath('/Users/yue_han/temp_for_mess/Plasmid_crosstalk_modeling/processed_data')
addpath('/Users/yue/Documents/Plasmid_crosstalk_modeling-/processed_data')
%% Further processing of mRNA degradation data 

% Load processed data 
gain_75_deg_data_table = readtable('degradation_curve_gain_75.xlsx','Sheet','single_sample_deg');
gain_75_deg_data_timeVec = readtable('degradation_curve_gain_75.xlsx','Sheet','timeVec'); 
gain_75_init_cond = [10,50,100,500,1000,5000,10000];

gain_70_deg_data_table = readtable('degradation_curve_gain_70.xlsx','Sheet','single_sample_deg');
gain_70_deg_data_timeVec = readtable('degradation_curve_gain_70.xlsx','Sheet','timeVec'); 
gain_70_init_cond = [10,50,100,500,1000,5000,10000,20000];

% Impute weird data using two approaches 
% plot_flag = true; 
% [imputed_data_gain70,sim_imputed_data_gain70,imputed_init_cond_gain70,sim_imputed_init_cond_gain70,full_sim_timeVec_gain70] = ...
%     extrap_mRNA_deg_cont(gain_70_deg_data_timeVec,gain_70_deg_data_table,gain_70_init_cond,plot_flag);
% [imputed_data_gain75,sim_imputed_data_gain75,imputed_init_cond_gain75,sim_imputed_init_cond_gain75,full_sim_timeVec_gain75] = ...
%     extrap_mRNA_deg_cont(gain_75_deg_data_timeVec,gain_75_deg_data_table,gain_75_init_cond,plot_flag);

% % Save imputed processed data 
% save('/Users/yue/Documents/Plasmid_crosstalk_modeling-/processed_data/imputed_RNAdeg_data_v2.mat'); 

%% Use the imputed actual time zero data for calibration curve 
load('/Users/yue/Documents/Plasmid_crosstalk_modeling-/processed_data/imputed_RNAdeg_data_v2.mat'); 

    % Perform linear regression
gain_70_lr_v2 = [sim_imputed_init_cond_gain70' \ sim_imputed_data_gain70(1,~isnan(sim_imputed_data_gain70(1,:)))',0];
gain_75_lr_v2 = [sim_imputed_init_cond_gain75' \ sim_imputed_data_gain75(1,~isnan(sim_imputed_data_gain75(1,:)))',0];

%     % Save to existing calibration curve as v2 
% save('/Users/yue/Documents/Plasmid_crosstalk_modeling-/calibration_curve/standard_curve_gain_70.mat',"gain_70_lr_v2",'-append'); 
% save('/Users/yue/Documents/Plasmid_crosstalk_modeling-/calibration_curve/standard_curve_gain_75.mat',"gain_75_lr_v2",'-append'); 