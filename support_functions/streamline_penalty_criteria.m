function satisfy_flag_matrix = streamline_penalty_criteria(all_high_res_penalty_terms,penalty_cutoff)

    % Condense high resolution penalty terms into 9 criteria and generate a
    % matrix indicating whether each of the 9 is satisfied for each sample 

    %% INPUT
    %:all_high_res_penalty_terms: # samples * 32 matrix 
    %:penalty_cutoff: 1 * 9 vector defining threshold criteria

    %% OUTPUT 
    %:satisfy_flag_matrix: # samples * 9 logical matrix indicating whether
    % a certain criteria is met in a sample 

    %% function
    
    baseline_dev = all_high_res_penalty_terms(:,1); 
    crosstalk_ratio_dev = all_high_res_penalty_terms(:,2); 
    baseline_penalty = all_high_res_penalty_terms(:,3:6); 
    
    T7_strong_positive_crosstalk_penalty_idx_list = 7:4:15;
    T7_strong_negative_crosstalk_penalty_idx_list = 19:4:27;
    T7_weak_positive_crosstalk_penalty_idx_list = 8:4:16;
    T7_weak_negative_crosstalk_penalty_idx_list = 20:4:28;
    sigma70_strong_positive_crosstalk_penalty_idx_list = 9:4:17;
    sigma70_strong_negative_crosstalk_penalty_idx_list = 21:4:29;
    sigma70_weak_positive_crosstalk_penalty_idx_list = 10:4:18;
    sigma70_weak_negative_crosstalk_penalty_idx_list = 22:4:30;
    T7_strong_crosstalk_penalty = all_high_res_penalty_terms(:,[T7_strong_positive_crosstalk_penalty_idx_list,T7_strong_negative_crosstalk_penalty_idx_list]); 
    T7_weak_crosstalk_penalty = all_high_res_penalty_terms(:,[T7_weak_positive_crosstalk_penalty_idx_list,T7_weak_negative_crosstalk_penalty_idx_list]);     
    sigma70_strong_crosstalk_penalty = all_high_res_penalty_terms(:,[sigma70_strong_positive_crosstalk_penalty_idx_list,sigma70_strong_negative_crosstalk_penalty_idx_list]); 
    sigma70_weak_crosstalk_penalty = all_high_res_penalty_terms(:,[sigma70_weak_positive_crosstalk_penalty_idx_list,sigma70_weak_negative_crosstalk_penalty_idx_list]); 
    
    large_positive_crosstalk_penalty = all_high_res_penalty_terms(:,end-1);
    residual_mRNA_penalty = all_high_res_penalty_terms(:,end); 
    
    %     % Reorganize penalty cutoff
    T7_strong_crosstalk_penalty_cutoff = penalty_cutoff([T7_strong_positive_crosstalk_penalty_idx_list,T7_strong_negative_crosstalk_penalty_idx_list]); 
    T7_weak_crosstalk_penalty_cutoff = penalty_cutoff([T7_weak_positive_crosstalk_penalty_idx_list,T7_weak_negative_crosstalk_penalty_idx_list]); 
    sigma70_strong_crosstalk_penalty_cutoff = penalty_cutoff([sigma70_strong_positive_crosstalk_penalty_idx_list,sigma70_strong_negative_crosstalk_penalty_idx_list]); 
    sigma70_weak_crosstalk_penalty_cutoff = penalty_cutoff([sigma70_weak_positive_crosstalk_penalty_idx_list,sigma70_weak_negative_crosstalk_penalty_idx_list]); 

        % Generate satisfy flag matrix that is #samples * 9
    satisfy_flag_matrix = [baseline_dev <= penalty_cutoff(1),crosstalk_ratio_dev <= penalty_cutoff(2),all(baseline_penalty <= penalty_cutoff(3:6),2),...
        all(T7_strong_crosstalk_penalty <= T7_strong_crosstalk_penalty_cutoff,2),all(T7_weak_crosstalk_penalty <= T7_weak_crosstalk_penalty_cutoff,2),...
        all(sigma70_strong_crosstalk_penalty <= sigma70_strong_crosstalk_penalty_cutoff,2),all(sigma70_weak_crosstalk_penalty <= sigma70_weak_crosstalk_penalty_cutoff,2),...
        large_positive_crosstalk_penalty <= penalty_cutoff(end - 1),residual_mRNA_penalty <= penalty_cutoff(end)];







end