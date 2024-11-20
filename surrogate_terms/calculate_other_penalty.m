function penalty = calculate_other_penalty(simulated_data,species_name_list)
    main_species_list = {'kanR','empty','GFP'};
    penalty = 0; 
    for dataset_idx = 1:length(simulated_data)
        simulated_data_single = simulated_data{dataset_idx}; 
        % Find all the RNA species 
        for main_species_idx = 1:length(main_species_list)
            main_species = main_species_list{main_species_idx}; 
            RNA_species_idx_list = contains(species_name_list,main_species) & contains(species_name_list,'RNA ') & ~contains(species_name_list,'RNase'); 
            % selected_RNA_species_name_list = species_name_list(RNA_species_idx_list); 
            filtered_simulated_data = simulated_data_single(end,RNA_species_idx_list); % Find the end of the time course
            penalty = penalty + sum(filtered_simulated_data); 
        end
    end
end