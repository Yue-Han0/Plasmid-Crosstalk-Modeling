function run_fit_qualitative_characteristics_single_v2(obj_idx,iter)


    %% Set correct path to avoid confusion of config files
    currentpath = pwd;
    addpath(genpath(sprintf('%s',currentpath)));
    rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath));
    if ~exist('param_est_run_save','dir')
        mkdir('param_est_run_save')
    end

    group_description = {'PE_T7_strong_no_empty','PE_T7_strong_empty','PE_T7_strong_empty_T7','PE_T7_strong_empty_sigma70',...
        'PE_T7_weak_no_empty','PE_T7_weak_empty','PE_T7_weak_empty_T7','PE_T7_weak_empty_sigma70',...
        'PE_sigma70_strong_no_empty','PE_sigma70_strong_empty','PE_sigma70_strong_empty_T7','PE_sigma70_strong_empty_sigma70',...
        'PE_sigma70_weak_no_empty','PE_sigma70_weak_empty','PE_sigma70_weak_empty_T7','PE_sigma70_weak_empty_sigma70'};
    param_info_path = 'param_info/parameters_v2.xlsx';
    
    problemObject = setProblemObject_v3(group_description,[],param_info_path); 
    dosing_information = create_dosing_info_from_problemObject(problemObject); 
    mode = 'PE';
    num_conc = 7;
    num_promotor = 4; 
    num_params = length(problemObject.Estimated);
    param_bounds = [problemObject.Estimated.Bounds]; 
    param_bounds = reshape(param_bounds,[2,length(problemObject.Estimated)]);
    param_bounds = param_bounds'; 

    % Preassign tracking variables
    options = optimoptions('fmincon','Display','iter'); 

    rng(iter)
    % Get initial point
    init_params = nan(1,num_params); 
    for est_param_idx = 1:length(problemObject.Estimated)
        log_init_param =  problemObject.Estimated(est_param_idx).TransformedBounds(1) + ...
            (problemObject.Estimated(est_param_idx).TransformedBounds(2) -  problemObject.Estimated(est_param_idx).TransformedBounds(1)) * rand;
        init_param =  exp(log_init_param);
        init_params(est_param_idx) = init_param; 
    end

    [estimated_params,fval,exitflag,output] = ...
        fmincon(@(x) wrapper_calculate_obj_single(x,problemObject,dosing_information,obj_idx,mode,num_conc,num_promotor),...
        init_params,[],[],[],[],param_bounds(:,1)',param_bounds(:,2)',[],options);


    save(sprintf('param_est_run_save/20240823_fit_single_qual_obj%d_iter%d.mat',obj_idx,iter)); 


end


