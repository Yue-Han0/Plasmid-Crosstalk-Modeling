function converted = log_conversion(large_number)
    % Convert the large number into log form and avoid issues with negative
    % values 

    if large_number > 0 
        converted = log10(large_number);
    elseif large_number < 0
        converted = -1 * log10(-1 * large_number); 
    else
        converted = 0; 
    end


end