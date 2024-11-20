function table_for_fit = remove_neg_val(table_for_fit)
    %% Remove negative values from training data table 

    % Keep track of the rows to be kept 
    keep_idx_list = []; 
    table_for_fit_vars = table_for_fit.Properties.VariableNames; 
    if any(strcmp(table_for_fit_vars,'mRNA_concentration'))
        for i = 1:length(table_for_fit.mRNA_concentration)
            if table_for_fit.mRNA_concentration(i) >= 0 || table_for_fit.Time(i) < 50
                keep_idx_list = [keep_idx_list,i];
                % If time zero concentration is less than zero, replace
                % with zero instead of removing the row 
                if table_for_fit.Time(i) < 50 && table_for_fit.mRNA_concentration(i) < 0
                    table_for_fit.mRNA_concentration(i) = 0;   
                end
            end
            
        end
    elseif any(strcmp(table_for_fit_vars,'GFP_concentration'))
        for i = 1:length(table_for_fit.GFP_concentration)
            if table_for_fit.GFP_concentration(i) >= 0 || table_for_fit.Time(i) < 50
                keep_idx_list = [keep_idx_list,i];
                if table_for_fit.Time(i) < 50 && table_for_fit.GFP_concentration(i) < 0
                    table_for_fit.GFP_concentration(i) = 0;   
                end
            end
        end
    end
    table_for_fit = table_for_fit(keep_idx_list,:);


end