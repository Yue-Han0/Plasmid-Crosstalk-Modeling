function [Time,Data] = process_crosstalk_data_from_source(data,source,option)
    

    switch source
        % If data source is a grouped data object, from constructing
        % fitProblemStruct 
        case 'grouped_data'
            all_group_numbers = data.Group; 
            all_group_number = unique(all_group_numbers);

            % Preassign 
            Time = cell(length(all_group_number),1);
            Data = cell(length(all_group_number),1); 

            for i = 1:length(all_group_number)

                % Select data of a data group 
                group_number = all_group_number(i); 
                group_data = data(data.Group == group_number,:); 

                % Convert into table and extract timeVec & time-course 
                data_table = groupedData2table(group_data); 
                Time{i,1} = data_table.Time;
                try 
                    GFP_concentration = data_table.GFP_concentration;
                    if all(GFP_concentration == 0)
                        mRNA_flag = true;
                        GFP_flag = false;
                    end
                catch
                    mRNA_flag = true; 
                    GFP_flag = false;
                end
                try
                    mRNA_concentration = data_table.mRNA_concentration;
                    if all(mRNA_concentration == 0)
                        GFP_flag = true; 
                        mRNA_flag = false;
                    end
                catch
                    GFP_flag = true; 
                    mRNA_flag = false;
                end
                
                % Assign the non-zero time-course
                if mRNA_flag
                    Data{i,1} = mRNA_concentration;
                elseif GFP_flag
                    Data{i,1} = GFP_concentration;
                end
                Time{i,1} = data_table.Time; 

            end

        % If data source is simulated data, from running fitted(fitResult) 
        case 'simulated_data'
            
            % Preassign 
            Time = cell(length(data),1);
            Data = cell(length(data),1); 

            for i = 1:length(data)
                Time{i,1} = data(i).Time;
                DataNames = data(i).DataNames; 

                % Get time-course for GFP and broc
                GFP_idx = strcmp(DataNames,'protein sfGFP*'); 
                GFP_time_course = data(i).Data(:,GFP_idx);

                broc_idx = strcmp(DataNames,'RNA utrbroc--no_protein'); 
                broc_time_course = data(i).Data(:,broc_idx);

                % Decide between sfGFP vs. broc as reporter 
                if isempty(GFP_time_course) || all(abs(GFP_time_course - 0) < 1)
                    broc_flag = true; 
                    sfGFP_flag = false; 
                elseif isempty(broc_time_course) || all(abs(broc_time_course - 0) < 1)
                    broc_flag = false; 
                    sfGFP_flag = true; 
                end

                % Assign the non-zero time-course
                if ~exist('option','var') % Only ouptut reporter time-course unless otherwise specified 
                    if broc_flag
                        Data{i,1} = broc_time_course;
                    elseif sfGFP_flag
                        Data{i,1} = GFP_time_course;
                    end
                elseif strcmp(option,'keep_all')
                    Data{i,1} = data(i).Data; 
                end
            end

        % If data source is from running custom simFunction 
        case 'sim_function'
            
    
    end

end