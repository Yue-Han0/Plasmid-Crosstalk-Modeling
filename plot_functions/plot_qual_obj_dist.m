function plot_qual_obj_dist(high_res_penalty_terms,penalty_cutoff,save_path)

    % Visualize penalty terms in scatter plots

    %% INPUT
    % :high_res_penalty_terms: a # parameter sets * # penalty terms matrix
    % containing qualitative objective values 
    % :penalty_cutoff: 
    % :save_path: an optional string output defining save path

    % Break down into (1) baseline deviation (2) crosstalk ratio deviation
    % (3) baseline deviation (4) T7 strong crosstalk (5) T7 weak crosstalk
    % (6) sigma70 strong crosstalk (7) sigma70 weak crosstalk (8) larger
    % positive crosstalk (9) residual mRNA concentration 
    all_marker_type = {'o','^','square','diamond'}; 
    colormap_name = 'hsv'; 
    

    figure; 
    
    % (1) baseline deviation 
    subplot(3,3,1)
    baseline_dev = high_res_penalty_terms(:,1); 
    baseline_dev_cutoff = penalty_cutoff(1); 
    c = 1:length(baseline_dev) + 1; % Add one duplicat datapoint to prevent endpoint dumping 
    scatter(1:length(baseline_dev) + 1,[baseline_dev;baseline_dev(end)],[],c) 
    colormap(gca,colormap_name)
    yline(baseline_dev_cutoff,'LineWidth',2,'Color',"#EDB120")
    title('Baseline Deviation')

    % (2) crosstalk ratio deviation 
    subplot(3,3,2)
    crosstalk_ratio_dev = high_res_penalty_terms(:,2); 
    crosstalk_ratio_dev_cutoff = penalty_cutoff(2); 
    scatter(1:length(crosstalk_ratio_dev) + 1,[crosstalk_ratio_dev;crosstalk_ratio_dev(end)],[],c)
    colormap(gca,colormap_name)
    yline(crosstalk_ratio_dev_cutoff,'LineWidth',2,'Color',"#EDB120")
    title('Crosstalk Ratio Deviation')

    % (3) Baseline penalty 
    subplot(3,3,3)
    baseline_penalty = high_res_penalty_terms(:,3:6); 
    baseline_penalty_cutoff = penalty_cutoff(3:6); 
    for prom_idx = 1:length(baseline_penalty_cutoff)
        scatter(1:size(baseline_penalty,1) + 1,[baseline_penalty(:,prom_idx);baseline_penalty(end,prom_idx)],[],c,all_marker_type{prom_idx})
        hold on
    end
    colormap(gca,colormap_name)
    yline(baseline_penalty_cutoff(prom_idx),'LineWidth',2,'Color',"#EDB120")
    title('Baseline Penalty Deviation')
    legend('T7 strong','T7 weak','sigma70 strong','sigma70 weak','Penalty Cutoff')

    % (4) T7 strong crosstalk 
    subplot(3,3,4)
    T7_strong_positive_crosstalk_penalty_idx_list = 7:4:15;
    T7_strong_negative_crosstalk_penalty_idx_list = 19:4:27;
    T7_strong_crosstalk_penalty = high_res_penalty_terms(:,[T7_strong_positive_crosstalk_penalty_idx_list,T7_strong_negative_crosstalk_penalty_idx_list]); 
    T7_strong_crosstalk_penalty_cutoff = penalty_cutoff([T7_strong_positive_crosstalk_penalty_idx_list,T7_strong_negative_crosstalk_penalty_idx_list]); 
    for case_idx = 1:3
        scatter(1:size(T7_strong_crosstalk_penalty,1) + 1,[T7_strong_crosstalk_penalty(:,case_idx);T7_strong_crosstalk_penalty(end,case_idx)],[],c,all_marker_type{case_idx}) % Positive crosstalk 
        hold on
        scatter(1:size(T7_strong_crosstalk_penalty,1) + 1,[T7_strong_crosstalk_penalty(:,3 + case_idx);T7_strong_crosstalk_penalty(end,3 + case_idx)],[],c,all_marker_type{case_idx},'filled') % Negative crosstalk
    end
    colormap(gca,colormap_name)
    yline(T7_strong_crosstalk_penalty_cutoff(prom_idx),'LineWidth',2,'Color',"#EDB120")
    title('T7 Strong Crosstalk')
    legend('Empty - Positive Crosstalk','Empty - Negative Crosstalk','Empty T7 - Positive Crosstalk','Empty T7 - Negative Crosstalk','Empty sigma70 - Positive Crosstalk',...
        'Empty sigma70 - Negative Crosstalk','Penalty Cutoff')

    % (5) T7 weak crosstalk 
    subplot(3,3,5)
    T7_weak_positive_crosstalk_penalty_idx_list = 8:4:16;
    T7_weak_negative_crosstalk_penalty_idx_list = 20:4:28;
    T7_weak_crosstalk_penalty = high_res_penalty_terms(:,[T7_weak_positive_crosstalk_penalty_idx_list,T7_weak_negative_crosstalk_penalty_idx_list]); 
    T7_weak_crosstalk_penalty_cutoff = penalty_cutoff([T7_weak_positive_crosstalk_penalty_idx_list,T7_weak_negative_crosstalk_penalty_idx_list]); 
    for case_idx = 1:3
        scatter(1:size(T7_weak_crosstalk_penalty,1) + 1,[T7_weak_crosstalk_penalty(:,case_idx);T7_weak_crosstalk_penalty(end,case_idx)],[],c,all_marker_type{case_idx}) % Positive crosstalk 
        hold on
        scatter(1:size(T7_weak_crosstalk_penalty,1) + 1,[T7_weak_crosstalk_penalty(:,3 + case_idx);T7_weak_crosstalk_penalty(end,3 + case_idx)],[],c,all_marker_type{case_idx},'filled') % Negative crosstalk
    end
    colormap(gca,colormap_name)
    yline(T7_weak_crosstalk_penalty_cutoff(prom_idx),'LineWidth',2,'Color',"#EDB120")
    title('T7 weak Crosstalk')
    legend('Empty - Positive Crosstalk','Empty - Negative Crosstalk','Empty T7 - Positive Crosstalk','Empty T7 - Negative Crosstalk','Empty sigma70 - Positive Crosstalk',...
        'Empty sigma70 - Negative Crosstalk','Penalty Cutoff')

    % (6) sigma70 strong crosstalk 
    subplot(3,3,6)
    sigma70_strong_positive_crosstalk_penalty_idx_list = 9:4:17;
    sigma70_strong_negative_crosstalk_penalty_idx_list = 21:4:29;
    sigma70_strong_crosstalk_penalty = high_res_penalty_terms(:,[sigma70_strong_positive_crosstalk_penalty_idx_list,sigma70_strong_negative_crosstalk_penalty_idx_list]); 
    sigma70_strong_crosstalk_penalty_cutoff = penalty_cutoff([sigma70_strong_positive_crosstalk_penalty_idx_list,sigma70_strong_negative_crosstalk_penalty_idx_list]); 
    for case_idx = 1:3
        scatter(1:size(sigma70_strong_crosstalk_penalty,1) + 1,[sigma70_strong_crosstalk_penalty(:,case_idx);sigma70_strong_crosstalk_penalty(end,case_idx)],[],c,all_marker_type{case_idx}) % Positive crosstalk 
        hold on
        scatter(1:size(sigma70_strong_crosstalk_penalty,1) + 1,[sigma70_strong_crosstalk_penalty(:,3 + case_idx);sigma70_strong_crosstalk_penalty(end,3 + case_idx)],[],c,all_marker_type{case_idx},'filled') % Negative crosstalk
    end
    colormap(gca,colormap_name)
    yline(sigma70_strong_crosstalk_penalty_cutoff(prom_idx),'LineWidth',2,'Color',"#EDB120")
    title('sigma70 strong Crosstalk')
    legend('Empty - Positive Crosstalk','Empty - Negative Crosstalk','Empty T7 - Positive Crosstalk','Empty T7 - Negative Crosstalk','Empty sigma70 - Positive Crosstalk',...
        'Empty sigma70 - Negative Crosstalk','Penalty Cutoff')

    % (7) sigma70 weak crosstalk 
    subplot(3,3,7)
    sigma70_weak_positive_crosstalk_penalty_idx_list = 10:4:18;
    sigma70_weak_negative_crosstalk_penalty_idx_list = 22:4:30;
    sigma70_weak_crosstalk_penalty = high_res_penalty_terms(:,[sigma70_weak_positive_crosstalk_penalty_idx_list,sigma70_weak_negative_crosstalk_penalty_idx_list]); 
    sigma70_weak_crosstalk_penalty_cutoff = penalty_cutoff([sigma70_weak_positive_crosstalk_penalty_idx_list,sigma70_weak_negative_crosstalk_penalty_idx_list]); 
    for case_idx = 1:3
        scatter(1:size(sigma70_weak_crosstalk_penalty,1) + 1,[sigma70_weak_crosstalk_penalty(:,case_idx);sigma70_weak_crosstalk_penalty(end,case_idx)],[],c,all_marker_type{case_idx}) % Positive crosstalk 
        hold on
        scatter(1:size(sigma70_weak_crosstalk_penalty,1) + 1,[sigma70_weak_crosstalk_penalty(:,3 + case_idx);sigma70_weak_crosstalk_penalty(end,3 + case_idx)],[],c,all_marker_type{case_idx},'filled') % Negative crosstalk
    end
    colormap(gca,colormap_name)
    yline(sigma70_weak_crosstalk_penalty_cutoff(prom_idx),'LineWidth',2,'Color',"#EDB120")
    title('sigma70 weak Crosstalk')
    legend('Empty - Positive Crosstalk','Empty - Negative Crosstalk','Empty T7 - Positive Crosstalk','Empty T7 - Negative Crosstalk','Empty sigma70 - Positive Crosstalk',...
        'Empty sigma70 - Negative Crosstalk','Penalty Cutoff')

    % (8) larger positive crosstalk 
    subplot(3,3,8)
    larger_positive_crosstalk_penalty = high_res_penalty_terms(:,end-1); 
    larger_positive_crosstalk_penalty_cutoff = penalty_cutoff(end-1); 
    scatter(1:length(larger_positive_crosstalk_penalty) + 1,[larger_positive_crosstalk_penalty;larger_positive_crosstalk_penalty(end)],[],c)
    colormap(gca,colormap_name)
    yline(larger_positive_crosstalk_penalty_cutoff,'LineWidth',2,'Color',"#EDB120")
    title('Larger Positive Crosstalk Penalty')

    % (9) residual mRNA concentration 
    subplot(3,3,9)
    residual_mRNA_penalty = high_res_penalty_terms(:,end); 
    residual_mRNA_penalty_cutoff = penalty_cutoff(end); 
    scatter(1:length(residual_mRNA_penalty) + 1,[residual_mRNA_penalty;residual_mRNA_penalty(end)],[],c)
    colormap(gca,colormap_name)
    yline(residual_mRNA_penalty_cutoff,'LineWidth',2,'Color',"#EDB120")
    title('Residual mRNA Penalty')


    if exist('save_path','var')
        save(gcf,save_path)
    end



end