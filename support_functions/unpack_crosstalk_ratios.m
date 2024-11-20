function varargout = unpack_crosstalk_ratios(crosstalk_ratio,promotor_name_list)
    

    vararg_var_names = cell(size(promotor_name_list)); 
    for prom_idx = 1:length(promotor_name_list)
        variable_name = sprintf('%s_crosstalk_ratio',promotor_name_list{prom_idx}); 
        prom_crosstalk_ratio = crosstalk_ratio(:,prom_idx); 
        % Concatenate prom_crosstalk_ratios from 3*1 cell to 3 * #dose
        % matrix 
        prom_crosstalk_ratio_vec = cell2mat(prom_crosstalk_ratio); 
        prom_crosstalk_ratio_mat = reshape(prom_crosstalk_ratio_vec,[length(prom_crosstalk_ratio_vec) / size(crosstalk_ratio,1),size(crosstalk_ratio,1)]); 
        eval(sprintf('varargout{prom_idx} = prom_crosstalk_ratio_mat;')); 
        vararg_var_names{prom_idx} = variable_name; 

    end




end