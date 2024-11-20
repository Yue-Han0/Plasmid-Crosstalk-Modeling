function [timeVec,species_time_course] = get_species_time_course(simData,species_oi)
    % Get the timeVec and time-course data for specific species

    %%% Input %%%
    % simData: a Simbiology simData object
    % species_oi: a (1 * nSpecies) vector of species names in simData

    %%% Output %%%
    % timeVec: a (nT * 1) vector
    % species_time_course: a (nT * nSpecies) vector

    species_idx_list = nan(length(species_oi),1); 
    timeVec = simData.Time; 
    species_time_course = nan(length(timeVec),length(species_oi)); 
    for idx = 1:length(species_oi)
        species_name = species_oi{idx};
        species_idx = find(strcmp(simData.DataNames,species_name)); 
        if isempty(species_idx)
            warning('%s not found in simData',species_name);
        else
            species_idx_list(idx) = species_idx; 
            species_time_course(:,idx) = simData.Data(:,species_idx); 
        end
    end

end