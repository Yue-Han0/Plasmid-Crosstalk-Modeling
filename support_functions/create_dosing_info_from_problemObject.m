function dosing_information_for_simFun = create_dosing_info_from_problemObject(problemObject)
    % Create the dosing information required for simFunction 

    % Construct a dosing information matrix based on problemObject.Doses
    dosing_information_mat = zeros(size(problemObject.Doses));
    for dose_idx = 1:size(problemObject.Doses,1)
        for target_idx = 1:size(problemObject.Doses,2)
            single_dose_amount = problemObject.Doses(dose_idx,target_idx).Amount;
            if ~isempty(single_dose_amount)
                dosing_information_mat(dose_idx,target_idx) = single_dose_amount; 
            end
        end
    end

    % Convert into a table for simFun
    dosing_information = cell(size(dosing_information_mat)); 
    for row_idx = 1:size(dosing_information_mat,1)
        for col_idx = 1:size(dosing_information_mat,2)
            dosing_table = array2table([0,dosing_information_mat(row_idx,col_idx)]); 
            dosing_table.Properties.VariableNames = {'Time','Amount'};
            dosing_information{row_idx,col_idx} = dosing_table; 
        end
    end

    dosing_information_for_simFun = dosing_information; 

end