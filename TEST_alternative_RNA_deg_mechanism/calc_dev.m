%% Define an error function 
function dev = calc_dev(params,simFunction,dosing_information,tStart,tEnd,RNA_degradation_data_table)

    [simTime,simData] = simFunction(params,tEnd,dosing_information,tStart:tEnd);
    dev = 0; 
    for conc_idx = 1:length(simData)
        simulated_timeVec = simTime{conc_idx}(2:end); % Remove first value 
        simulated_data = simData{conc_idx}(2:end,1);
        exp_timeVec = RNA_degradation_data_table.Time(RNA_degradation_data_table.Group == conc_idx);
        exp_data = RNA_degradation_data_table.mRNA_concentration(RNA_degradation_data_table.Group == conc_idx);

        % Remove duplicate values from simulated data 
        [simulated_timeVec_unique,unique_idx,~] = unique(simulated_timeVec); 
        simulated_data_unique = simulated_data(unique_idx);

        % Extrapolate data 
        simulated_data_extrap = interp1(simulated_timeVec_unique,simulated_data_unique,exp_timeVec,'linear','extrap');

        % Calculate deviation 
        dev = dev + sum(sum((simulated_data_extrap - exp_data).^2)); 
    end
end