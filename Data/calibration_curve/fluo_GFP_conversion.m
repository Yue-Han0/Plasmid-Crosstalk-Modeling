function converted = fluo_GFP_conversion(fluo,GFP,gain,slope,intercept)
    % Conversion between fluorescence and GFP at a specific gain  
        % fluo - a matrix of fluorescence value 
        % GFP - a matrix of GFP concentration values
        % gain - a value representing plate reader gain 
        % slope, intercept - values that can be input to override the
        % loaded numbers 

    if ~isempty(fluo) && ~isempty(GFP)
        error('Can only do one-way conversion!')
    elseif isempty(fluo) && isempty(GFP)
        error('Missing data input')
    end
    
    if ~exist('slope','var') || ~exist('intercept','var')
        if isequal(gain,50)
            standard_curve_file = load('calibration_curve/sfGFP_standard_curve_gain_50.mat');
            lr_coeff = standard_curve_file.gain_50_lr; 
        elseif isequal(gain,75)
            standard_curve_file = load('calibration_curve/sfGFP_standard_curve_gain_75.mat');
            lr_coeff = standard_curve_file.gain_75_lr;
        end
        slope = lr_coeff(1);
        intercept = lr_coeff(2); 
    end
    if isempty(fluo) % converting GFP to fluorescence 
        converted = slope .* GFP + intercept; 
    elseif isempty(GFP) % converting fluorescence to GFP 
        converted = (fluo - intercept) ./ slope; 
    end


end