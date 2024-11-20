clear
clc

currentpath = pwd; 
addpath(genpath(currentpath))
addpath('/Users/yue/Documents/Plasmid_crosstalk_modeling-/processed_data/')
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath))

num_iter = 120; 

%% Build model 
Model_struct = sbioloadproject('simbiology_models/RNA_deg.sbproj');
Model = Model_struct.m1; 

% Let's first add new species
RNase_bound = addspecies(Model,'RNase_bound');
mRNA_deact = addspecies(Model,'mRNA_deact');

% (1) For each RNase binding reaction, change the product into RNase_bound
% (2) Duplicate the RNase degradation reaction and add deactivated mRNA as a product 
reactions_to_be_added = {}; 
for rxn_idx = 1:length(Model.Reactions)

    reaction_selected = Model.Reactions(rxn_idx); 
    chemical_reaction = reaction_selected.Reaction; 
    divider_idx = strfind(chemical_reaction,'->'); 
    rxn_lhs = chemical_reaction(1:divider_idx); 
    rxn_rhs = chemical_reaction(divider_idx + 1:end); 
    % find RNase binding reactions and replace product 
    if contains(rxn_lhs,' RNase ') && contains(rxn_rhs,':RNase')  

        % Get product name and remove product 
        product_name = get(reaction_selected.product); 
        rmproduct(reaction_selected, product_name.Name); 

        % Replace with RNase_bound
        addproduct(reaction_selected,RNase_bound,1); 
    end

    % Find RNA degradation reaction and (1) replace reactant (2) get ready to duplicate those 
    % (3) Remove AGMP/CUMP 
    if contains(rxn_lhs,':RNase') && contains(rxn_rhs,' RNase ')

        reactions_to_be_added{end + 1} = reaction_selected;

        % Get reactant name and remove reactant 
        reactant_name = get(reaction_selected.Reactants);
        rmreactant(reaction_selected,reactant_name.Name); 

        % Replace with RNase_bound
        addreactant(reaction_selected,RNase_bound,1); 

        % Remove AGMP/CUMP 
        rmproduct(reaction_selected,'AGMP');
        rmproduct(reaction_selected,'CUMP'); 

    end

end

% Add a new reaction to account for RNase binding to deactivated mRNA 
    % Add paramet
TXTL_RNAdeg_R_proc = addparameter(Model,'TXTL_RNAdeg_kc_proc',1); 

  % Add new reactions
RNase_binding_mRNAdeact = addreaction(Model, 'mRNA_deact + RNase <-> RNase_bound');
RNase_binding_mRNAdeact_kineticlaw = addkineticlaw(RNase_binding_mRNAdeact,'MassAction');
RNase_binding_mRNAdeact_kineticlaw.ParameterVariableNames = {'TXTL_RNAdeg_F','TXTL_RNAdeg_R'};

% Duplciate the reactions and add deactivated mRNA as a product 
for dup_rxn_idx = 1:length(reactions_to_be_added)
    rxn_oi = reactions_to_be_added{dup_rxn_idx};
    dup_rxn_oi = copyobj(rxn_oi,Model); 

    % Add product to the duplicated reaction 
    addproduct(dup_rxn_oi,mRNA_deact,1); 

end

% Use different kinetic parameters for processive and non-processive RNA
% degradation 
    % These should act on the same set of reactions 
for proc_deg_rxn_idx = 1:length(reactions_to_be_added)
    rxn_oi = reactions_to_be_added{proc_deg_rxn_idx};

    % Change the associated parameter variable names 
    rxn_oi.KineticLaw.ParameterVariableNames = {'TXTL_RNAdeg_kc_proc'}; 
end

% Clean up - remove duplicate reactions, remove unused species
rxn_idx = 1;
while rxn_idx < length(Model.Reactions)
    rxn_oi = Model.Reactions(rxn_idx); 
    rxn_name = rxn_oi.Reaction; 
    all_model_rxn_name = {Model.Reactions.Reaction}; 
    rxn_name_idx_list = find(strcmp(rxn_name,all_model_rxn_name)); 
    if length(rxn_name_idx_list) > 1 
        rxn_to_be_deleted = Model.Reactions(rxn_name_idx_list(2)); 
        delete(rxn_to_be_deleted); 
    else
        rxn_idx = rxn_idx + 1; 
    end
end

   % Remove species that are not used in any reactions
reaction_name_list = {Model.Reactions.Reaction};
species_idx = 1; 
while species_idx <= length(Model.Species)
    selected_species = Model.Species(species_idx); 
    exist_flag = false; 
    for rxn_idx = 1:length(reaction_name_list)
        rxn_name = reaction_name_list{rxn_idx}; 
        if contains(rxn_name,selected_species.Name)
            exist_flag = true; 
        end
    end
    if ~exist_flag
        delete(selected_species)
    else
        species_idx = species_idx + 1; 
    end
end
    
    % For newly added parameters/reactions, add parameters in Reaction.KineticLaw to model parameters and delete
% them from kinetic law object 
Mobj_rxns = get(Model,'Reactions');
for rxn_idx = 1:length(Mobj_rxns)
    rxn = Mobj_rxns(rxn_idx);
    try
        rxn_params = rxn.KineticLaw.Parameters;
    catch
        rxn_params = []; 
    end
    for rxn_param_idx = 1:length(rxn_params)
        rxn_param = rxn_params(rxn_param_idx);
        rxn_param_name = rxn_param.Name; 
    
        if ~strcmp(rxn_param_name,rxn_param.Name)
            % Change parameter variable names
            parameter_name_idx = strcmp(rxn.KineticLaw.ParameterVariableNames,rxn_param.Name); 
            rxn.KineticLaw.ParameterVariableNames{parameter_name_idx} = rxn_param_name;
        end
        % Is this parameter already in model parameters?
        global_param = sbioselect(Model.Parameters,'Name',rxn_param_name);
        if ~isempty(global_param)
            % if already in model parameters, select and delete
            fprintf('\n %s already exists in model object,deleting',rxn_param_name)
            delete(rxn_param); 
        else
            % If doesn't exist in model parameters, move to model
            % parameters 
                % Make sure parameter name matches 
            if ~strcmp(rxn_param_name,rxn_param.Name)
                %If doesn't match, set parameter name
                set(rxn_param,'Name',rxn_param_name); 
            end
            test_obj = move(rxn_param,Model);
        end
    end
end

%% Load experimental data 
    % From processed data directly 
% gain_70_timeVec = readtable('processed_data/degradation_curve_gain_70.xlsx','Sheet','timeVec');
% gain_70_RNA_time_course = readtable('processed_data/degradation_curve_gain_70.xlsx','Sheet','degradation');
% gain_75_timeVec = readtable('processed_data/degradation_curve_gain_75.xlsx','Sheet','timeVec');
% gain_75_RNA_time_course = readtable('processed_data/degradation_curve_gain_75.xlsx','Sheet','degradation');
% 
%         % Get initial concentration 
% RNA_init_conc_gain70 = gain_70_RNA_time_course{1,2:9}; 
% RNA_init_conc_gain75 = gain_75_RNA_time_course{1,2:8}; 
%         % Get timeVec and time-course 
% RNA_timeVec_gain70 = gain_70_timeVec{:,2}; 
% RNA_timeVec_gain75 = gain_75_timeVec{:,2}; 
% RNA_fluo_time_course_gain70 = gain_70_RNA_time_course{2:end,11:18}; 
% RNA_fluo_time_course_gain75 = gain_75_RNA_time_course{2:end,10:16}; 

    % From imputed processed data 
imputed_processed_RNA_deg_file = load('imputed_RNAdeg_data.mat'); 
RNA_fluo_time_course_gain70 = imputed_processed_RNA_deg_file.imputed_data_gain70; 
RNA_init_conc_gain70 = imputed_processed_RNA_deg_file.imputed_init_cond_gain70; 
RNA_timeVec_gain70 = imputed_processed_RNA_deg_file.gain_70_deg_data_timeVec{:,2}; 
% RNA_fluo_time_course_gain70 = imputed_processed_RNA_deg_file.sim_imputed_data_gain70; 
% RNA_init_conc_gain70_eff  = imputed_processed_RNA_deg_file.sim_imputed_init_cond_gain70;

        % Convert fluorescence to mRNA concentration based on calibration
        % curve
% Note: these two should technically give the same mRNA concentrations, but
% due to experimental noise and fitting erros there's a small difference
% (~10%) between them. Accounting for uncertainty propogation here would be
% insane, so let's use gain 70 to fit degradation. 
RNA_fluo_time_course_gain70 = RNA_fluo_time_course_gain70(:,~any(isnan(RNA_fluo_time_course_gain70),1));
RNA_time_course_gain70 = fluo_mRNA_conversion(RNA_fluo_time_course_gain70,[],70); 

    % Add a time delay to timeVec and data
% RNA_timeVec_gain70 = RNA_timeVec_gain70 + 4 * 60; % 4-min time delay 
% RNA_timeVec_gain70 = [0;RNA_timeVec_gain70]; % Add a true time zero 
% RNA_time_course_gain70 = [RNA_init_conc_gain70;RNA_time_course_gain70]; % Add init conc for the true time zero 

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


%% Build a fitProblem struct 
RNA_deg_probObject = fitproblem(); % Initialize a fitProblem object
        % Data
RNA_deg_probObject.Data = groupedData(RNA_degradation_data_table); % define a groupedData object 
RNA_deg_probObject.Data.Properties.IndependentVariableName = 'Time'; 
RNA_deg_probObject.Data.Properties.GroupVariableName = 'Group'; 
        % Model 
RNA_deg_probObject.Model = Model; 
RNA_deg_probObject.ResponseMap = '[RNA utrbroc--no_protein] = mRNA_concentration';
        % Doses
RNA_deg_probObject.Doses = createDoses(RNA_deg_probObject.Data,{'mRNA_concentration_add'}); 
for dose_idx = 1:size(RNA_deg_probObject.Doses,1)
    RNA_deg_probObject.Doses(dose_idx,1).TargetName = 'RNA utrbroc--no_protein';
end
        % Estimates 
param_info_table = readtable('parameters_RNAdeg_newMech.xlsx');
Mobj_parameters = get(Model,'Parameters'); 
paramsToEstimate = {}; 
initialValues = []; 
bounds = []; 
for param_idx = 2:length(Mobj_parameters)
    Mobj_param = Mobj_parameters(param_idx);
    table_idx = find(strcmp(param_info_table.Name,Mobj_param.Name));
    paramsToEstimate{end + 1} = strcat('log(',Mobj_param.Name,')');
    initVal = table2array(param_info_table(table_idx,'InitVal')); 
    initialValues = [initialValues initVal];
    LB = table2array(param_info_table(table_idx,'LB')); 
    UB = table2array(param_info_table(table_idx,'UB')); 
    bounds = [bounds;[LB,UB]]; 
end
RNA_deg_probObject.Estimated = estimatedInfo(paramsToEstimate,'InitialValue',initialValues,'Bounds',bounds);
% RNA_deg_probObject.FunctionName = 'lsqcurvefit'; 
RNA_deg_probObject.Pooled = true; 
RNA_deg_probObject.ProgressPlot = 0; 

% [fitResults,~] = fit(RNA_deg_probObject); % A fit problem object 

% figure;
% for conc_idx = 1:5%size(RNA_time_course_gain70,2)
%     plot(RNA_timeVec_gain70,RNA_time_course_gain70(:,conc_idx))
%     hold on 
% end

% % Run fit 
% all_fitResults = cell(num_iter,1); 
% parfor iter = 1:num_iter
%     % Modify initial conditions, sample within bounds
%     rng(iter)
%     temp_problemObject = RNA_deg_probObject; 
%     % Use initial parameter values for first run 
%     % if ~isequal(iter,1)
%         % Then sample initial parameter values within parameter bounds
%     for est_param_idx = 1:length(temp_problemObject.Estimated)
%         temp_problemObject.Estimated(est_param_idx).InitialTransformedValue =  temp_problemObject.Estimated(est_param_idx).TransformedBounds(1) +...
%             (temp_problemObject.Estimated(est_param_idx).TransformedBounds(2) -  temp_problemObject.Estimated(est_param_idx).TransformedBounds(1)) * rand;
%         temp_problemObject.Estimated(est_param_idx).InitialValue = exp(temp_problemObject.Estimated(est_param_idx).InitialTransformedValue); 
%     end
%     % end   
%     try
%         [fitResults,~] = fit(temp_problemObject); % A fit problem object 
%         all_fitResults{iter,1} = fitResults;
%     catch
%         fprintf('Fitting #%d failed',iter)
%     end
% 
% end
% 
% all_SSE = nan(num_iter,1); 
% for iter = 1:num_iter
%     if ~isempty(all_fitResults{iter})
%         fitResult = all_fitResults{iter,1}; 
%         all_SSE(iter) = fitResult.SSE;
%     end
% end
% 
% % Check fitting 
% [~,opt_idx] = min(all_SSE); 
% opt_fitResult = all_fitResults{opt_idx}; 
% opt_fitted_data = fitted(opt_fitResult); 
% for conc_idx = 1:length(opt_fitted_data)
%     simData_single = opt_fitted_data(conc_idx);
%     subplot(ceil(sqrt(length(opt_fitted_data))),ceil(sqrt(length(opt_fitted_data))),conc_idx)
%     plot(simData_single.Time(2:end),simData_single.Data(2:end,1),'LineWidth',1.5,'Color','r')
%     hold on 
%     plot(RNA_degradation_data_table.Time(RNA_degradation_data_table.Group == conc_idx),...
%         RNA_degradation_data_table.mRNA_concentration(RNA_degradation_data_table.Group == conc_idx),'LineWidth',1.5,'Color','k')
% end

%% Build a simFunction and use sbiofit 
    % Create simFunction 
model_species = {Model.Species.Name}; 
species_oi_name = 'RNA utrbroc--no_protein'; 
species_oi_idx = find(strcmp(species_oi_name,model_species)); 
track_species = [model_species(1:species_oi_idx - 1),model_species(species_oi_idx + 1:end)];
simFunction = create_simFun_from_problemObject(RNA_deg_probObject,track_species); 
dosing_information= create_dosing_info_from_problemObject(RNA_deg_probObject); 
tStart = 0; 
tEnd = 2940; 
options = optimoptions('fmincon','Display','iter');

    % Preassign 
all_fval = nan(num_iter,1);
all_estimated_params = nan(num_iter,length(RNA_deg_probObject.Estimated)); 

    % Run Parameter Estimation 
for iter = 1:num_iter
    log_transformed_init_params = nan(1,length(RNA_deg_probObject.Estimated));
    init_params = nan(1,length(RNA_deg_probObject.Estimated));
    for est_param_idx = 1:length(RNA_deg_probObject.Estimated)
        log_transformed_init_params(est_param_idx) =  RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(1) +...
            (RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(2) -  RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(1)) * rand;
        init_params(est_param_idx) = exp(log_transformed_init_params(est_param_idx)); 
    end
    [estimated_params,fval,exitflag,output] = fmincon(@(x)calc_dev(x,simFunction,dosing_information,tStart,tEnd,RNA_degradation_data_table),init_params,[],[],[],[],bounds(:,1),bounds(:,2),[],options);

    all_fval(iter) = fval; 
    all_estimated_params(iter,:) = estimated_params; 
end

    % Plot out SSE distirbution 
figure;
histogram(log10(all_fval))
xlabel('log_{10}SSE')
ylabel('# Occurence')

    % Plot out simulated vs. experimental
[~,sort_idx] = sort(all_fval,'ascend'); 
opt_estimated_params = all_estimated_params(sort_idx(1),:); 
[simTime,simData] = simFunction(opt_estimated_params,tEnd,dosing_information,tStart:tEnd); 
figure; 
for conc_idx = 1:length(simData)
    subplot(ceil(sqrt(length(simData))),ceil(sqrt(length(simData))),conc_idx)
    plot(simTime{conc_idx},simData{conc_idx}(:,1),'LineWidth',1.5,'Color','r')
    hold on 
    plot(RNA_degradation_data_table.Time(RNA_degradation_data_table.Group == conc_idx),...
        RNA_degradation_data_table.mRNA_concentration(RNA_degradation_data_table.Group == conc_idx),'LineWidth',1.5,'Color','k')
end

    % Check parameter distribution to assess parameter sloppiness 
figure; 
lower_SSE_idx_list = (all_fval < 1e+08); 
lower_SSE_params = all_estimated_params(lower_SSE_idx_list,:); 
higher_SSE_params = all_estimated_params(~lower_SSE_idx_list,:); 
param_names = {RNA_deg_probObject.Estimated.Name};
for param_idx = 1:length(RNA_deg_probObject.Estimated)
    
    subplot(ceil(sqrt(length(RNA_deg_probObject.Estimated))),ceil(sqrt(length(RNA_deg_probObject.Estimated))),param_idx)
    histogram(log(lower_SSE_params(:,param_idx)))
    hold on
    histogram(log(higher_SSE_params(:,param_idx)))
    title(strrep(param_names{param_idx},'_',' '))

end

    % Define an error function 
function dev = calc_dev(params,simFunction,dosing_information,tStart,tEnd,RNA_degradation_data_table)

    [simTime,simData] = simFunction(params,tEnd,dosing_information,tStart:tEnd);
    dev = 0; 
    for conc_idx = 1:length(simData)
        simulated_timeVec = simTime{conc_idx}(2:end); % Remove first value 
        simulated_data = simData{conc_idx}(2:end,1);
        exp_timeVec = RNA_degradation_data_table.Time(RNA_degradation_data_table.Group == conc_idx);
        exp_data = RNA_degradation_data_table.mRNA_concentration(RNA_degradation_data_table.Group == conc_idx);

        % Remove duplicate values from simulated data 
        [simulated_timeVec_unique,unique_idx,~] = unique(simulated_timeVec); 
        simulated_data_unique = simulated_data(unique_idx);

        % Extrapolate data 
        simulated_data_extrap = interp1(simulated_timeVec_unique,simulated_data_unique,exp_timeVec,'linear','extrap');

        % Calculate deviation 
        dev = dev + sum(sum((simulated_data_extrap - exp_data).^2)); 
    end
end
% [simulated_time_ctrl,simulated_data_ctrl] = simFunction_ctrl(estimated_params,tEnd,dosing_information,tStart:tEnd);

    % sbiofit 
% ResponseMap = '[RNA utrbroc--no_protein] = mRNA_concentration';
% for iter = 1:48
%     estimated_info = RNA_deg_probObject.Estimated;
%     for est_param_idx = 1:length(estimated_info)
%         estimated_info(est_param_idx).InitialTransformedValue =  RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(1) +...
%             (RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(2) -  RNA_deg_probObject.Estimated(est_param_idx).TransformedBounds(1)) * rand;
%         estimated_info(est_param_idx).InitialValue = exp(RNA_deg_probObject.Estimated(est_param_idx).InitialTransformedValue); 
%     end
%     try
%         sbiofit_fitResults = sbiofit(Model,RNA_deg_probObject.Data,ResponseMap,estimated_info,RNA_deg_probObject.Doses); 
%     catch
%         fprintf('Run # %d failed',iter)
%     end
% end


%% Check if the mRNA calibration curve is actually similar in alternative mechanisms 

%     % By manually writing ODE and fit slope and intercept
%     % simultnaeously
%     % Load mRNA degrdation data
% gain_50_timeVec_file= readcell('processed_data/degradation_curve_gain_50.xlsx',sheet='timeVec'); 
% gain_50_data_file = readcell('processed_data/degradation_curve_gain_50.xlsx',sheet='degradation');
% 
% gain_100_timeVec_file= readcell('processed_data/degradation_curve_gain_100.xlsx',sheet='timeVec'); 
% gain_100_data_file = readcell('processed_data/degradation_curve_gain_100.xlsx',sheet='degradation');
% 
% gain_70_timeVec_file= readcell('processed_data/degradation_curve_gain_70.xlsx',sheet='timeVec'); 
% gain_70_data_file = readcell('processed_data/degradation_curve_gain_70.xlsx',sheet='degradation');
% 
% gain_75_timeVec_file= readcell('processed_data/degradation_curve_gain_75.xlsx',sheet='timeVec'); 
% gain_75_data_file = readcell('processed_data/degradation_curve_gain_75.xlsx',sheet='degradation');
% 
%     % Convert to matrix forms ready for fitting 
% gain_50_timeVec = cell2mat(gain_50_timeVec_file(2:end,2)); 
% gain_100_timeVec = cell2mat(gain_100_timeVec_file(2:end,2)); 
% gain_70_timeVec = cell2mat(gain_70_timeVec_file(2:end,2)); 
% gain_75_timeVec = cell2mat(gain_75_timeVec_file(2:end,2)); 
%         % get initial mRNA concentrations (Gain 50 & 100) 
% mRNA_concentration_str_vec = gain_100_data_file(1,2:end); % this should be the same for gain 50 and gain 100
% mRNA_concentration_vec = nan(1,length(mRNA_concentration_str_vec)); 
% for conc_idx = 1:length(mRNA_concentration_vec)
%     conc_str = mRNA_concentration_str_vec{conc_idx}; 
%     conc_num = conc_str(1:end-3); 
%     mRNA_concentration_vec(1,conc_idx) = str2double(conc_num);  
% end
%     % get initial mRNA concentrations (Gain 70 & 75) 
% mRNA_concentration_70_updated_vec = cell2mat(gain_70_data_file(1,2:9));
% mRNA_concentration_75_updated_vec = cell2mat(gain_75_data_file(1,2:8));
% 
%     % gain data 
% gain_50_data = cell2mat(gain_50_data_file(2:end,2:end)); 
% gain_100_data = cell2mat(gain_100_data_file(2:end,2:end)); 
% 
% gain_70_data = cell2mat(gain_70_data_file(2:end,11:end));
% gain_75_data = cell2mat(gain_75_data_file(2:end,10:end));
% 
% time_delay = 240; 
%     % Shift timeVec to 4 minutes back (240s) 
% actual_gain50_timeVec = gain_50_timeVec + time_delay; 
% actual_gain100_timeVec = gain_100_timeVec + time_delay; 
% 
% actual_gain70_timeVec = gain_70_timeVec + time_delay;
% actual_gain75_timeVec = gain_75_timeVec + time_delay; 
%     % Add an actual time zero at the start 
% actual_gain50_timeVec = [0;actual_gain50_timeVec];
% actual_gain100_timeVec = [0;actual_gain100_timeVec];
% 
% actual_gain70_timeVec = [0;actual_gain70_timeVec];
% actual_gain75_timeVec = [0;actual_gain75_timeVec];
%     % choose a degradation kinetics to fit to 
% %         % 1st-order degradation 
% % kinetics = 3; 
% %         % used 1st pass curve as initial point & 0.5x and 2x as lb and ub 
% % init_params = [4,1e-06,0.001]; % Assume first 2 variables are slope and intercept, the 3rd is k_deg in 1st-order kinetics 
% % init_params_lb = [0,0,0.00001];
% % init_params_ub = [16,1,1]; 
% 
% % init_params = [17.45,0,7e-01,1000];
% % init_params_lb = [8,0,1e-02,1000];
% % init_params_ub = [32,10,10,1e07];
% 
%     % Mixed RNA degradation mechanism
% init_params = [17.45,0,1e04,1e07,1e-04,1e02,1e-04]; 
%         % params: (3) k_binding (4) k_unbinding (5) k_deg (6) RNase0 (7)
%         % k_deg_proc
% init_params_lb = [0,-1e-08,1e03,1e06,1e-05,1e02,1e-05];
% init_params_ub = [16,1e-08,1e06,1e10,1e-01,1e05,1e-01];
%     % Define variables 
% % experimental_data = struct('timeVec',actual_gain100_timeVec,'data',gain_100_data,'real_init',mRNA_concentration_vec);
% % gain = 100; 
% % fmincon_options = optimoptions('fmincon','Display','iter-detailed');
% experimental_data = struct('timeVec',actual_gain70_timeVec,'data',gain_70_data,'real_init',mRNA_concentration_70_updated_vec);
% gain = 70; 
% fmincon_options = optimoptions('fmincon','Display','iter-detailed');
% 
% estimated_params_all = cell(48,1);
% fval_all = nan(48,1); 
% 
% for iter = 1:480
%     if isequal(iter,1)
%         iter_init_params = init_params; 
%     else
%         iter_init_params = init_params_lb + rand(size(init_params)) .* (init_params_ub - init_params_lb); 
%     end
%     [estimated_params,fval] = fmincon(@(x) calc_dev(x,experimental_data,gain),iter_init_params,[],[],[],[],init_params_lb,init_params_ub,[],fmincon_options); 
%     estimated_params_all{iter,1} = estimated_params;
%     fval_all(iter) = fval; 
% end
% 
% [~,best_idx] = min(fval_all); 
% best_estimated_params = estimated_params_all{best_idx,1}; 
% [sorted_fval_all,sort_idx] = sort(fval_all,'ascend'); 
% sorted_estimated_params_all = estimated_params_all(sort_idx,:);
% 
% % Plot out fitted results & experimental data
% tspan = [experimental_data.timeVec(1),experimental_data.timeVec(end)];
% figure; 
% for init_conc_idx = 1:length(experimental_data.real_init)
% 
%     subplot(ceil(sqrt(length(experimental_data.real_init))),ceil(sqrt(length(experimental_data.real_init))),init_conc_idx)
% 
%     init_cond = [experimental_data.real_init(init_conc_idx),best_estimated_params(7),0,0]; 
%     ode_options = odeset('NonNegative',1:length(init_cond));
%     [ode_timeVec_opt,ode_sol_opt] = ode15s(@(t,y)build_mRNA_deg_ode_v2(t,y,best_estimated_params),tspan,init_cond,ode_options);
%     raw_predicted_data = ode_sol_opt(:,1); 
%     fluo_data = fluo_mRNA_conversion([],raw_predicted_data,70,best_estimated_params(1),best_estimated_params(2));
% 
%     plot(ode_timeVec_opt,fluo_data,'LineWidth',1.5,'Color','r')
%     hold on 
%     plot(experimental_data.timeVec(2:end),experimental_data.data(:,init_conc_idx),'LineWidth',1.5,'Color','k')
% 
% end
% 
% function error = calc_dev(params,experimental_data,gain)
% 
%     % experimental_data: a data structure with timeVec, data, and real_init
%         % timeVec: nT * 1 matrix with time in s 
%         % real_init: 1 * num_init_conc matrix representing added mRNA concentration 
%         % data: nT * num_init_conc matrix representing time-course
%         % degradation of mRNA 
%     init_fluo = fluo_mRNA_conversion([],experimental_data.real_init,gain,params(1),params(2)); 
%     if ~all(isequal(init_fluo,experimental_data.data(1,:)))
%         experimental_data.data = [init_fluo;experimental_data.data];
%     end
% 
%     error = 0; 
%     for init_conc_idx = 1:length(experimental_data.real_init)
% 
% 
%         % Filter the experimental dataset, discard timepoints after a
%         % negative value occurs 
%         tspan = [experimental_data.timeVec(1),experimental_data.timeVec(end)]; % ode solves for the entire time period
% 
%         negative_idx_list = find(experimental_data.data(:,init_conc_idx)<=0);
%         if ~isempty(negative_idx_list)
%             valid_timeVec = experimental_data.timeVec(1:negative_idx_list(1)-1);
%             valid_fluorescence_data = experimental_data.data(1:negative_idx_list(1)-1,init_conc_idx); 
%         else
%             valid_timeVec = experimental_data.timeVec;
%             valid_fluorescence_data = experimental_data.data(:,init_conc_idx); 
%         end
% 
%         % if initial mRNA concentration is larger than 500nM,
%         % fluorescence keeps increasing after time 0, so assume that the
%         % peak represents mRNA initial concentrations 
%         if experimental_data.real_init(init_conc_idx) >= 500 
%             valid_timeVec = valid_timeVec(2:end);
%             valid_fluorescence_data = valid_fluorescence_data(2:end); 
%         end
% 
%         % Assign initial conditions for different reaction mechanisms 
%         init_cond = [experimental_data.real_init(init_conc_idx),params(7),0,0]; 
%         ode_options = odeset('NonNegative',1:length(init_cond));
% 
%         % Solve for mRNA degradation curve 
%         [ode_timeVec,ode_sol] = ode15s(@(t,y)build_mRNA_deg_ode_v2(t,y,params),tspan,init_cond,ode_options);
%         raw_predicted_data = ode_sol(:,1); 
% 
%         % Extrapolate/intrapolate to valid timeVec for comparison with
%         % experimental data 
%         if isequal(ode_timeVec(end),tspan(end))
%             predicted_data = interp1(ode_timeVec,raw_predicted_data,valid_timeVec,'linear','extrap'); 
%             predicted_data_fluo = fluo_mRNA_conversion([],predicted_data,gain,params(1),params(2)); 
%             error = error + sum((valid_fluorescence_data - predicted_data_fluo).^2); 
%         else
%             error = error + 1e08;
%         end
%     end
% end
% 
% function dydt = build_mRNA_deg_ode_v2(~,y,params) % Build ODE for mixed processive and non-processive degradation 
%     % y(1) - functional mRNA; y(2) - Free RNase; y(3) - Bound RNase; y(4) -
%     % deactivated mRNA 
% 
%     % params - (1) slope (2) Intercept (3) k_binding (4) k_unbinding (5) k_deg_nonProc (6)
%     % k_deg_proc (7) RNase_0 
% 
%     dydt = nan(4,1);
% 
%     % Reaction rates 
%     v_RNase_binding_mRNAf = params(3) * y(1) * y(2);
%     v_RNase_unbinding_mRNAf = params(4) * y(3); 
%     v_RNase_binding_mRNAde = params(3) * y(4) * y(2); 
%     v_RNase_unbinding_mRNAde = params(4) * y(3);
%     v_proc_deg = params(6) * y(3);
%     v_nonProc_deg = params(5) * y(3);
% 
%     % Mass balance
%     dydt(1) = v_RNase_unbinding_mRNAf - v_RNase_binding_mRNAf;
%     dydt(2) = v_RNase_unbinding_mRNAf + v_RNase_unbinding_mRNAde - v_RNase_binding_mRNAf - v_RNase_binding_mRNAde + v_proc_deg + v_nonProc_deg; 
%     dydt(3) = v_RNase_binding_mRNAf + v_RNase_binding_mRNAde - v_RNase_unbinding_mRNAf - v_RNase_unbinding_mRNAde  - v_proc_deg - v_nonProc_deg; 
%     dydt(4) = v_RNase_unbinding_mRNAde - v_RNase_binding_mRNAde + v_nonProc_deg;
% 
% end
    