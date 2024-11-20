function gen_feasibility_constraints_value()
    % Generate the goal of each qualitative objective terms 

    % For high-resolution penalty 
    sampling_result_file = load('test_save_files/202408_high_res_penalty_valid_terms.mat','valid_high_res_penalty_terms'); 
    high_res_penalty_terms = sampling_result_file.valid_high_res_penalty_terms; 
    
    % Qunatitative values such as baseline and crosstalk ratio deviation
    % and residue RNA are constrained using median values 
    median_baseline_penalty_value = median(high_res_penalty_terms(:,1),'omitmissing');
    median_crosstalk_ratio_penalty_value = median(high_res_penalty_terms(:,2),'omitmissing');
    median_other_penalty_value = median(high_res_penalty_terms(:,end),'omitmissing'); 
    
    % Positive Crosstalk and Negative crosstalk are constrained with a 0.05
    % to account for small deviation from 1 
    positive_crosstalk_goals = -0.05 .* ones(1,12);
    negative_crosstalk_goals = -0.05 .* ones(1,12); 
        % With the exception that T7 strong + empty sigma70 (#9) & T7 weak +
        % empty T7 (#6) do not need to comply to positive crosstalk 
    positive_crosstalk_goals(6) = 100;
    positive_crosstalk_goals(9) = 100; % Set these to be an arbitrary large value so it's always satisfied 


    penalty_goals = [median_baseline_penalty_value,median_crosstalk_ratio_penalty_value,...
        -1e-04,-1e-04,-1e-04,-1e-04,...
        positive_crosstalk_goals,negative_crosstalk_goals,-1e-04,median_other_penalty_value];

    save('test_save_files/high_resolution_penalty_goals.mat','penalty_goals'); 






end