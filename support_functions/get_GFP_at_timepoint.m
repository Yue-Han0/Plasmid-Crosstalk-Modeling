    
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
            GFP_val = data(min_idx,1);
        else
            % If not, interpolate to get the desired timepoint using a 90s
            % interval as the experimental data 
            interp_timeVec = 0:90:timeVec(end); 
            GFP_time_course_interp = interp1(timeVec,data(:,1),interp_timeVec,"linear"); 
            GFP_val = GFP_time_course_interp(interp_timeVec==timepoint_oi); 
    
        end
    elseif strcmp(timepoint_oi,'end')
        GFP_val = data(end,1); 
    end
    if isempty(GFP_val)
        disp('pause')
    end

end