function RNA_deg_probObject = build_RNAdeg_fitProblem_struct(kinetics,Model,Data)

    RNA_deg_probObject = fitproblem(); % Initialize a fitProblem object

        % Data
    RNA_deg_probObject.Data = groupedData(Data); % define a groupedData object 
    RNA_deg_probObject.Data.Properties.IndependentVariableName = 'Time'; 
    RNA_deg_probObject.Data.Properties.GroupVariableName = 'Group'; 
    
        % Doses
    RNA_deg_probObject.Doses = createDoses(RNA_deg_probObject.Data,{'mRNA_concentration_add'}); 
    for dose_idx = 1:size(RNA_deg_probObject.Doses,1)
        RNA_deg_probObject.Doses(dose_idx,1).TargetName = 'RNA utrbroc--no_protein';
    end

        % Model 
    RNA_deg_probObject.Model = Model; 
    RNA_deg_probObject.ResponseMap = '[RNA utrbroc--no_protein] = mRNA_concentration';

    RNA_deg_probObject.Pooled = true; 
    RNA_deg_probObject.ProgressPlot = 0; 

    switch kinetics
        case 'firstOrder'

            paramsToEstimate = {'log(k_deg)'}; 
            initialValues = 0.001; 
            bounds = [0.00001,1];
            RNA_deg_probObject.Estimated = estimatedInfo(paramsToEstimate,'InitialValue',initialValues,'Bounds',bounds);

        case 'MichaelisMenten'

            paramsToEstimate = {'log(Vm)','log(Km)'}; 
            initialValues = [7e-01,1000]; 
            bounds = [1e-02,1000;1e+03,1e+07;];
            RNA_deg_probObject.Estimated = estimatedInfo(paramsToEstimate,'InitialValue',initialValues,'Bounds',bounds);

        case 'MassAction'

                % Estimates 
            param_info_table = readtable('parameters_RNAdeg.xlsx');
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

        case {'MixedDeg','TestPrototype','TestBindingSites'}

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

    end
    
end