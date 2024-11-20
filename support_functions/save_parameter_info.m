function save_parameter_info(Mobj,output_path)


    parameter_names = cell(length(Mobj.Parameters),1);
    parameter_init_values = nan(length(Mobj.Parameters),1); 
    parameter_lb = nan(length(Mobj.Parameters),1); 
    parameter_ub = nan(length(Mobj.Parameters),1); 

    for i = 1:length(Mobj.Parameters)
        parameter_names{i,1} = Mobj.Parameters(i).Name; 
        parameter_init_values(i,1) = Mobj.Parameters(i).Value; 
        % Default lb to be 0.1-fold and ub 10-fold
        if contains(Mobj.Parameters(i).Name,'_F')
            % fix the binding to be 1 
                % !! NOTE: if _F becomes present in names this needs to be
                % modified 
            parameter_lb(i,1) = Mobj.Parameters(i).Value; 
            parameter_ub(i,1) = Mobj.Parameters(i).Value; 
        else
            parameter_lb(i,1) = Mobj.Parameters(i).Value * 0.1; 
            parameter_ub(i,1) = Mobj.Parameters(i).Value * 10; 
        end
    end
    
    parameter_table = table(parameter_names,parameter_init_values,parameter_lb,parameter_ub,...
        'VariableNames',{'Name','InitVal','LB','UB'});
    writetable(parameter_table,output_path); 
end