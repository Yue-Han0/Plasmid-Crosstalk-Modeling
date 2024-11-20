function all_penalty_terms = wrapper_calcQualObj_from_params(params,problemObject,dosing_information)
    % Given a set of parameters for the plasmid crosstalk master model,
    % solve the time-course and calculate the qualitative objective terms 

        %%  Create problemObject 
            
    %     % For now let's work on protein-level crosstalk only 
    % group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
    %     'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
    %     'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
    %     'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
    % param_info_path = 'parameters_unbounded.xlsx'; % Use unbounded parameters 
    % problemObject = setProblemObject_v3(group_description,[],param_info_path); 
    % dosing_information = create_dosing_info_from_problemObject(problemObject); 
        
    
        %% Evaluate penalty values  
    [all_penalty_terms,~,~,~] = wrapper_calculate_obj_higher_resolution(params,...
        problemObject,dosing_information,[],[],'PE',7,4);


end