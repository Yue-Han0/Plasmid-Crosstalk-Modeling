function converted = fluo_mRNA_conversion(fluo,mRNA,gain,slope,intercept)
    % Conversion between fluorescence and mRNA at a specific gain  
        % fluo - a matrix of fluorescence value 
        % mRNA - a matrix of mRNA values
        % gain - a value representing plate reader gain 
        % slope, intercept - values that can be input to override the
        % loaded numbers 

    if ~isempty(fluo) && ~isempty(mRNA)
        error('Can only do one-way conversion!')
    elseif isempty(fluo) && isempty(mRNA)
        error('Missing data input')
    end
    
    if ~exist('slope','var') || ~exist('intercept','var')
        if isequal(gain,70)
            standard_curve_file = load('calibration_curve/standard_curve_gain_70.mat');
            lr_coeff = standard_curve_file.gain_70_lr_v2;
        elseif isequal(gain,75)
            standard_curve_file = load('calibration_curve/standard_curve_gain_75.mat');
            lr_coeff = standard_curve_file.gain_75_lr_v2;
        end
        slope = lr_coeff(1);
        intercept = lr_coeff(2); 
    end
    if isempty(fluo) % converting mRNA to fluorescence 
        converted = slope .* mRNA + intercept; 
    elseif isempty(mRNA) % converting fluorescence to mRNA 
        converted = (fluo - intercept) ./ slope; 
    end


end