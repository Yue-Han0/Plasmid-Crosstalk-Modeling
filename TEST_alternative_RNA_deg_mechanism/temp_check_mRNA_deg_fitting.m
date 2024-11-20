% Check mRNA production and degradation analysis 
resultFileName1 = 'param_est_run_save/20231120_param_est_run1127_mRNA_deg.mat'; 
resultFileName2 = 'param_est_run_save/20231127_param_est_run1450_mRNA_deg.mat'; 
resultFile1 = load(resultFileName1); 
resultFile2 = load(resultFileName2); 

    % Calculate Km, Vmax, and Vmax/Km in each case 
num_iter = length(resultFile1.all_fitResults); 

% For result file 1 
all_Km_resultFile1 = nan(num_iter,1); 
all_Vmax_resultFile1 = nan(num_iter,1); 
all_Vmax_over_Km_resultFile1 = nan(num_iter,1); 
all_SSE_resultFile1 = nan(num_iter,1); 
all_estimated_params_1 = nan(num_iter,3); 

% For result file 2 
all_Km_resultFile2 = nan(num_iter,1); 
all_Vmax_resultFile2 = nan(num_iter,1); 
all_Vmax_over_Km_resultFile2 = nan(num_iter,1); 
all_SSE_resultFile2 = nan(num_iter,1);
all_estimated_params_2 = nan(num_iter,3); 

for iter = 1:num_iter
    fitResult_1 = resultFile1.all_fitResults{iter}; 
    for row_idx = 1:height(fitResult_1.ParameterEstimates)
        eval(sprintf('%s = %d;',fitResult_1.ParameterEstimates.Name{row_idx},fitResult_1.ParameterEstimates.Estimate(row_idx))); 
    end

    all_Km_resultFile1(iter) = TXTL_RNAdeg_R + TXTL_RNAdeg_kc; 
    all_Vmax_resultFile1(iter) = TXTL_RNAdeg_kc * RNase_0; 
    all_Vmax_over_Km_resultFile1(iter) = TXTL_RNAdeg_kc * RNase_0 / (TXTL_RNAdeg_R + TXTL_RNAdeg_kc); 

    all_SSE_resultFile1(iter) = fitResult_1.SSE; 

    all_estimated_params_1(iter,:) = [TXTL_RNAdeg_R,TXTL_RNAdeg_kc,RNase_0]; 

end

for iter = 1:num_iter
    fitResult_2 = resultFile2.all_fitResults{iter}; 
    for row_idx = 1:height(fitResult_2.ParameterEstimates)
        eval(sprintf('%s = %d;',fitResult_2.ParameterEstimates.Name{row_idx},fitResult_2.ParameterEstimates.Estimate(row_idx))); 
    end

    all_Km_resultFile2(iter) = TXTL_RNAdeg_R + TXTL_RNAdeg_kc; 
    all_Vmax_resultFile2(iter) = TXTL_RNAdeg_kc * RNase_0; 
    all_Vmax_over_Km_resultFile2(iter) = TXTL_RNAdeg_kc * RNase_0 / (TXTL_RNAdeg_R + TXTL_RNAdeg_kc); 

    all_SSE_resultFile2(iter) = fitResult_2.SSE; 

    all_estimated_params_2(iter,:) = [TXTL_RNAdeg_R,TXTL_RNAdeg_kc,RNase_0]; 

end

combined_SSE = [all_SSE_resultFile1;all_SSE_resultFile2];
combined_Vmax = [all_Vmax_resultFile1;all_Vmax_resultFile2];
combined_Km = [all_Km_resultFile1;all_Km_resultFile2];
combined_Vmax_over_Km = [all_Vmax_over_Km_resultFile1;all_Vmax_over_Km_resultFile2];
combined_all_estimated_params = [all_estimated_params_1;all_estimated_params_2];

% Plot out the SSE distribution 
figure; 
histogram(log10(combined_SSE));
xlabel('log_{10}(SSE)')
ylabel('# Occurence')
set(gca,'FontSize',14)

% Select SSE in that biggest peak 
selected_idx_list = (log10(combined_SSE) < 8); 

selected_Vmax = combined_Vmax(selected_idx_list); 
selected_Km = combined_Km(selected_idx_list); 
selected_Vmax_over_Km = combined_Vmax_over_Km(selected_idx_list); 
selected_estimated_params = combined_all_estimated_params(selected_idx_list,:); 

% parameter distribution 
    % Vmax/Km 
figure;
histogram(selected_Vmax_over_Km)
xlabel('V_{max}/K_m')
ylabel('# Occurence')
set(gca,'FontSize',14)

    % Vmax
figure;
histogram(selected_Vmax)
xlabel('V_{max}')
ylabel('# Occurence')
set(gca,'FontSize',14)

    % Km 
figure;
histogram(selected_Km)
xlabel('K_m')
ylabel('# Occurence')
set(gca,'FontSize',14)

    % TXTL_RNAdeg_R_broc
figure;
histogram(selected_estimated_params(:,1))
xlabel('TXTL RNAdeg R Broc')
ylabel('# Occurence')
title('RNA_{deg,R} estimated value distribution')
set(gca,'FontSize',14)

    % TXTL_RNAdeg_kc
figure;
histogram(selected_estimated_params(:,2))
xlabel('TXTL RNAdeg kc Broc')
ylabel('# Occurence')
title('RNA_{deg,kc,broc} estimated value distribution')
set(gca,'FontSize',14)

    % RNase0
figure;
histogram(selected_estimated_params(:,3))
xlabel('RNase0')
ylabel('# Occurence')
title('RNase_0 estimated value distribution')
set(gca,'FontSize',14)


%% Calculate 95% confidence interval for the 4 TXTL parameters to use as stricter bounds 

    % From plot above, it's reasonable to assume the parameter distribution
    % to be normal 
param_bounds_CI = nan(size(selected_estimated_params,2),2); 
for est_param_idx = 1:size(selected_estimated_params,2)

    single_selected_estimated_params = selected_estimated_params(:,est_param_idx); 
    pd = fitdist(single_selected_estimated_params,'Normal'); 
    ci = paramci(pd,'Alpha',.05);

    param_bounds_CI(est_param_idx,:) = ci(:,1)'; 

end



