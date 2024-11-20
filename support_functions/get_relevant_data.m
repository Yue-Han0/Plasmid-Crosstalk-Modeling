function [group_number,selected_table] = get_relevant_data(data_table,data_description,group_description)
    % group_description - a cell array with description of data to be
    % selected
    % group_number - a vector of the group number for selected data
    % selected_table - a section of data table requested 

    group_number = []; 
    for group_idx = 1:length(group_description)
        group_description_s = group_description{group_idx};
        group_number = [group_number;find(strcmp(group_description_s,data_description))]; 
    end
    selected_idx = []; 
    for group_num = group_number'
        selected_idx = [selected_idx;find(data_table.Group==group_num)];
    end
    selected_table = data_table(selected_idx,:); 

end