function crosstalk_ratio = calculate_crosstalk_ratio_v2(Time,Data,num_conc,num_promotor,mode)

    % Input
        % Time: a [num_promotor * num_conc * 4] * 1 cell array with a time
        % vector 
        % Data: a [num_promotor * num_conc * 4] * 1 cell array with
        % mRNA_concentration ('TX' mode) or GFP_concentration ('PE' mode) 
        % num_conc: # reporter plasmid concentrations tested 
        % num_promotor: # reporter promotor tested 
        % mode: a string from {'TX','PE'} standing for transcription-level
        % or protein-level crosstalk 
    % Output
        % crosstalk_ratio: num_conc * num_promotor cell of crosstalk
        % ratios, each cell containing 3 values representing crosstalk
        % ratio of empty, empty T7, empty sigma70 
    if strcmp(mode,'TX') && num_promotor > 3
        error('Too many promotors for transcription-level crosstalk')
    end
    num_combination = 4; 
    timepoint_oi = 10800; 
    % timepoint_oi = 'end'; 

    crosstalk_ratio = cell(3,num_promotor); 
    for promotor_idx = 1:num_promotor
        no_empty_GFP_allconc = nan(num_conc,1); 
        empty_GFP_allconc = nan(num_conc,1); 
        empty_T7_GFP_allconc = nan(num_conc,1); 
        empty_sigma70_GFP_allconc = nan(num_conc,1); 
        for conc_idx = 1:num_conc

            no_empty_baseline_data_idx = (promotor_idx - 1) * num_combination * num_conc + conc_idx;
            empty_data_idx = (promotor_idx - 1) * num_combination * num_conc + num_conc + conc_idx;
            empty_T7_data_idx = (promotor_idx - 1) * num_combination * num_conc + 2 * num_conc + conc_idx;
            empty_sigma70_data_idx = (promotor_idx - 1) * num_combination * num_conc + 3 * num_conc + conc_idx;

            no_empty_baseline_Time = Time{no_empty_baseline_data_idx,1}; 
            empty_Time = Time{empty_data_idx,1}; 
            empty_T7_Time = Time{empty_T7_data_idx,1}; 
            empty_sigma70_Time = Time{empty_sigma70_data_idx,1}; 

            no_empty_baseline_Data = Data{no_empty_baseline_data_idx,1}; 
            empty_Data = Data{empty_data_idx,1}; 
            empty_T7_Data = Data{empty_T7_data_idx,1}; 
            empty_sigma70_Data = Data{empty_sigma70_data_idx,1}; 
            
            if strcmp(mode,'PE')
                no_empty_GFP_allconc(conc_idx,1) = get_GFP_at_timepoint(no_empty_baseline_Time,no_empty_baseline_Data,timepoint_oi); 
                empty_GFP_allconc(conc_idx,1) = get_GFP_at_timepoint(empty_Time,empty_Data,timepoint_oi); 
                empty_T7_GFP_allconc(conc_idx,1) = get_GFP_at_timepoint(empty_T7_Time,empty_T7_Data,timepoint_oi); 
                empty_sigma70_GFP_allconc(conc_idx,1) = get_GFP_at_timepoint(empty_sigma70_Time,empty_sigma70_Data,timepoint_oi); 
            elseif strcmp(mode,'TX')
                no_empty_broc_allconc(conc_idx,1) = get_broc_at_timepoint(no_empty_baseline_Time,no_empty_baseline_Data); 
                empty_broc_allconc(conc_idx,1) = get_broc_at_timepoint(empty_Time,empty_Data); 
                empty_T7_broc_allconc(conc_idx,1) = get_broc_at_timepoint(empty_T7_Time,empty_T7_Data); 
                empty_sigma70_broc_allconc(conc_idx,1) = get_broc_at_timepoint(empty_sigma70_Time,empty_sigma70_Data); 
            end


        end
        if strcmp(mode,'PE')
            crosstalk_ratio{1,promotor_idx} = empty_GFP_allconc ./ no_empty_GFP_allconc; 
            crosstalk_ratio{2,promotor_idx} = empty_T7_GFP_allconc ./ no_empty_GFP_allconc;     
            crosstalk_ratio{3,promotor_idx} = empty_sigma70_GFP_allconc ./ no_empty_GFP_allconc; 
        elseif strcmp(mode,'TX')
            crosstalk_ratio{1,promotor_idx} = empty_broc_allconc ./ no_empty_broc_allconc; 
            crosstalk_ratio{2,promotor_idx} = empty_T7_broc_allconc ./ no_empty_broc_allconc;     
            crosstalk_ratio{3,promotor_idx} = empty_sigma70_broc_allconc ./ no_empty_broc_allconc; 
        end

    end
    
    function GFP_val = get_GFP_at_timepoint(timeVec,data,timepoint_oi)
        % Get the simulated time-course for GFP 
        % Find the closest timepoint in Time to timepoint_oi
        if isnumeric(timepoint_oi)
            timepoint_oi_vec = repmat(timepoint_oi,[length(timeVec),1]); 
            timepoint_diff = abs(timeVec - timepoint_oi_vec); 
            [min_diff,min_idx] = min(timepoint_diff); 
            if isequal(min_diff,0)
                % if timepoint_oi is in the vector, output the corresponding GFP
                % value 
                GFP_val = data(min_idx);
            else
                % If not, interpolate to get the desired timepoint using a 90s
                % interval as the experimental data 
                interp_timeVec = 0:90:timeVec(end); 
                GFP_time_course_interp = interp1(timeVec,data,interp_timeVec,"linear"); 
                GFP_val = GFP_time_course_interp(interp_timeVec==timepoint_oi); 
        
            end
        elseif strcmp(timepoint_oi,'end')
            GFP_val = data(end); 
        end
        if isempty(GFP_val)
            disp('pause')
        end

    end

    function broc_val = get_broc_at_timepoint(~,data)
        % Get the simulated time-course for GFP 
        [~,max_idx] = max(data); 
        broc_val = data(max_idx);

    end
end
