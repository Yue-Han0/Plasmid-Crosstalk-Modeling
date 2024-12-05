function plot_simulated(Time,Data,option,exp_Time,exp_Data,species_option_struct)

% Centralized function to plot out simulated results in various options 
    % Input
        % Time & Data:  N*1 cell containing time vectors and simulated reporter concentration for individual datasets
        % from process_crosstalk_data_from_source
        % option: {'baseline','crosstalk_ratio','crosstalk','exhaustive','ConcSpecies','BundledSpecies'}
        % exp_Time & exp_Data (optional): N * 1 cell containing time vectors and experimental reporter concentration for individual datasets
        % from process_crosstalk_data_from_source

% Note: Currently applies to protein-level crosstalk only. 

    % Plot time-course data in the crosstalk format
    promotor_name_list = {'T7_strong','T7_weak','sigma70_strong','sigma70_weak'};
    conc_vec = [0.5,1,2.5,5,10,15,30];
    all_colors = {'g',[0.5,0.5,0.5],'r','b'};
    num_promotors = length(Data)/28; 
    cmap = hsv(256);
    colorIndices = round(linspace(1, size(cmap, 1) - 64,length(conc_vec)));
    sampledColors = cmap(colorIndices, :);

    figure; 
    switch option
        case 'baseline'
            for prom_idx = 1:num_promotors
                promotor_name = promotor_name_list{prom_idx}; 
                subplot(2,2,prom_idx)
                for conc_idx = 1:length(conc_vec)
                    data_idx = (prom_idx - 1) * 28 + conc_idx; 
                    simulated_time = Time{data_idx,1}; 
                    simulated_data = Data{data_idx,1}; 
                    plot(simulated_time ./ 3600,simulated_data(:,1),'LineWidth',1.5,'Color',sampledColors(conc_idx,:))
                    hold on 
                end
                xlabel('Time(hr)')
                ylabel('Concentration(nM)')
                if contains(promotor_name,'sigma')
                    promotor_name_for_plot = strcat('\',strrep(promotor_name,'_',' '));
                else
                    promotor_name_for_plot = strrep(promotor_name,'_',' ');
                end
                title(promotor_name_for_plot)
                legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')
                set(gca,'FontSize',14)
            end
        case 'baseline_short'
            for prom_idx = 1:num_promotors
                promotor_name = promotor_name_list{prom_idx}; 
                subplot(2,2,prom_idx)
                for conc_idx = 1:length(conc_vec)
                    data_idx = (prom_idx - 1) * 28 + conc_idx; 
                    simulated_time = Time{data_idx,1}; 
                    simulated_data = Data{data_idx,1}; 
                    selected_simulated_time = simulated_time(simulated_time <= 18000); 
                    selected_simulated_data = simulated_data(simulated_time <= 18000); 
                    plot(selected_simulated_time./3600,selected_simulated_data(:,1),'LineWidth',1.5,'Color',sampledColors(conc_idx,:))
                    hold on 
                end
                xlabel('Time(hr)')
                ylabel('Concentration(nM)')
                if contains(promotor_name,'sigma')
                    promotor_name_for_plot = strcat('\',strrep(promotor_name,'_',' '));
                else
                    promotor_name_for_plot = strrep(promotor_name,'_',' ');
                end
                title(promotor_name_for_plot)
                if isequal(prom_idx,num_promotors)
                    legend('0.5nM','1nM','2.5nM','5nM','10nM','15nM','30nM')
                end
                set(gca,'FontSize',14)
            end
        case 'crosstalk_ratio'
            crosstalk_ratios = calculate_crosstalk_ratio_v2(Time,Data,length(conc_vec),num_promotors,'PE'); 
            for prom_idx = 1:num_promotors
                subplot(2,2,prom_idx)
                promotor_name = promotor_name_list{prom_idx}; 
                prom_crosstalk_ratio = crosstalk_ratios(:,prom_idx); 
                scatter(conc_vec,prom_crosstalk_ratio{1},'MarkerFaceColor',all_colors{2},'LineWidth',0.75)
                hold on 
                scatter(conc_vec,prom_crosstalk_ratio{2},'MarkerFaceColor',all_colors{3},'LineWidth',0.6)
                scatter(conc_vec,prom_crosstalk_ratio{3},'MarkerFaceColor',all_colors{4})
                plot(conc_vec,prom_crosstalk_ratio{1},'Color',all_colors{2},'LineWidth',1.7)
                hold on 
                plot(conc_vec,prom_crosstalk_ratio{2},'Color',all_colors{3},'LineWidth',1.6)
                plot(conc_vec,prom_crosstalk_ratio{3},'Color',all_colors{4},'LineWidth',1.5)
                if contains(promotor_name,'sigma')
                    promotor_name_for_plot = strcat('\',strrep(promotor_name,'_',' '));
                else
                    promotor_name_for_plot = strrep(promotor_name,'_',' ');
                end
                title(promotor_name_for_plot)
                yline(1,'LineWidth',1.5)
                if isequal(prom_idx,num_promotors)
                    legend('empty','empty T7','empty \sigma70','No Crosstalk')
                end
                xlabel('Plasmid Concentration (nM)')
                ylabel('Crosstalk Ratio')
                set(gca,'FontSize',14)
            end

        case 'exhaustive'
            for conc_idx = 1:length(conc_vec)
                plasmid_conc = conc_vec(conc_idx); 
                for prom_idx = 1:num_promotors
                    
                    plot_idx = (conc_idx - 1) * num_promotors + prom_idx;
                    subplot(length(conc_vec),num_promotors,plot_idx)
        
                    promotor_name = promotor_name_list{prom_idx};
        
                    if isequal(conc_idx,1)
                        title(strrep(promotor_name,'_',' '))
                        hold on 
                    end
                    if isequal(conc_idx,7)
                        xlabel('time(s)')
                        hold on 
                    end
                    if isequal(prom_idx,1)
                        ylabel(sprintf('%.1f nM',plasmid_conc))
                        hold on 
                    end
                       
                    for empty_idx = 1:4
                        data_idx = (prom_idx - 1) * 4 * length(conc_vec) + (empty_idx - 1) * length(conc_vec) + conc_idx; 
                        if exist('exp_Time','var') && exist('exp_Data','var')
                            experimental_time = exp_Time{data_idx}; 
                            experimental_data = exp_Data{data_idx}; 
                        end
                        simulated_time = Time{data_idx}; 
                        simulated_data = Data{data_idx}; 
                
                        plot(simulated_time,simulated_data(:,1),'LineWidth',1.5,'LineStyle','--','Color',all_colors{empty_idx})
                        hold on
                        if exist('exp_Time','var') && exist('exp_Data','var')
                            plot(experimental_time,experimental_data,'LineWidth',1.5,'LineStyle','-','Color',all_colors{empty_idx})
                            hold on 
                        end
        
                    end
                end
            end
        
            line1 = plot(nan,nan,'Color','g','LineStyle','-','DisplayName','no empty experimental');
            line2 = plot(nan,nan,'Color','g','LineStyle','--','DisplayName','no empty simulated');
            line3 = plot(nan,nan,'Color',[0.5,0.5,0.5],'LineStyle','-','DisplayName','empty experimental');
            line4 = plot(nan,nan,'Color',[0.5,0.5,0.5],'LineStyle','--','DisplayName','empty simulated');
            line5 = plot(nan,nan,'Color','b','LineStyle','-','DisplayName','empty sigma70 experimental');
            line6 = plot(nan,nan,'Color','b','LineStyle','--','DisplayName','empty sigma70 simulated');
            line7 = plot(nan,nan,'Color','r','LineStyle','-','DisplayName','empty T7 experimental');
            line8 = plot(nan,nan,'Color','r','LineStyle','--','DisplayName','empty T7 simulated');
            if exist('exp_Time','var') && exist('exp_Data','var')
                legend([line1,line2,line3,line4,line5,line6,line7,line8])
            else
                legend([line2,line4,line6,line8])
            end

        case 'ConcSpecies' % Plot out all species time-course for low concentration crosstalk - need to specify promotor of interest 
            
            % Sanity check: species option struct exists 
            if ~exist('species_option_struct','var')
                error('Specify dateset to be plotted')
            end

            % Sanity check: Promotor name and species name list are specified 
            if ~isfield(species_option_struct,'promotor_oi')
                error('No promotor specified!')
            end
            if ~isfield(species_option_struct,'species_name_list')
                error('Species Name List not specified')
            end
            if ~isfield(species_option_struct,'conc_idx')
                error('Dataset index not specified')
            end

            promotor_oi = species_option_struct.promotor_oi;
            species_name_list = species_option_struct.species_name_list; 
            conc_idx = species_option_struct.conc_idx; 

            % Sanity check: simulated data contains all species 
            num_species = length(species_name_list);
            sample_simData = Data{1};
            if ~isequal(size(sample_simData,2),num_species)
                error('Mismatch in number of species in model')
            end

            prom_idx = find(strcmp(promotor_name_list,promotor_oi)); 
                % No empty
            simulated_time_no_empty = Time{(prom_idx - 1) * 4 * length(conc_vec) + conc_idx}; 
            simulated_data_no_empty = Data{(prom_idx - 1) * 4 * length(conc_vec) + conc_idx}; 
                % Empty 
            simulated_time_empty = Time{(prom_idx - 1) * 4 * length(conc_vec) + 1 * length(conc_vec) + conc_idx}; 
            simulated_data_empty = Data{(prom_idx - 1) * 4 * length(conc_vec) + 1 * length(conc_vec) + conc_idx}; 
                % Empty T7 
            simulated_time_emptyT7 = Time{(prom_idx - 1) * 4 * length(conc_vec) + 2 * length(conc_vec) + conc_idx}; 
            simulated_data_emptyT7 = Data{(prom_idx - 1) * 4 * length(conc_vec) + 2 * length(conc_vec) + conc_idx}; 
                % Empty sigma 70 
            simulated_time_empty_sigma70 = Time{(prom_idx - 1) * 4 * length(conc_vec) + 3 * length(conc_vec) + conc_idx}; 
            simulated_data_empty_sigma70 = Data{(prom_idx - 1) * 4 * length(conc_vec) + 3 * length(conc_vec) + conc_idx}; 

            for species_idx = 1:num_species
                subplot(ceil(sqrt(num_species)),ceil(sqrt(num_species)),species_idx)

                % Plot out 4 cases 
                plot(simulated_time_no_empty,simulated_data_no_empty(:,species_idx),'LineWidth',1.5,'Color',all_colors{1})
                hold on
                plot(simulated_time_empty,simulated_data_empty(:,species_idx),'LineWidth',1.5,'Color',all_colors{2})
                plot(simulated_time_emptyT7,simulated_data_emptyT7(:,species_idx),'LineWidth',1.5,'LineStyle','--','Color',all_colors{4})
                plot(simulated_time_empty_sigma70,simulated_data_empty_sigma70(:,species_idx),'LineWidth',1.5,'LineStyle','-.','Color',all_colors{3})

                species_name = species_name_list{species_idx}; 
                species_name_for_plot = strrep(species_name,'_',' ');
                title(species_name_for_plot)

            end
            
        case 'BundledSpecies' % Bundle bound species for more intuitive visualization; also remove species that are close to zero and report 
           

            % Same sanity check as ConcSpecies case 
                % Sanity check: species option struct exists 
                if ~exist('species_option_struct','var')
                    error('Specify dateset to be plotted')
                end
    
                % Sanity check: Promotor name and species name list are specified 
                if ~isfield(species_option_struct,'promotor_oi')
                    error('No promotor specified!')
                end
                if ~isfield(species_option_struct,'species_name_list')
                    error('Species Name List not specified')
                end
                if ~isfield(species_option_struct,'conc_idx')
                    error('Dataset index not specified')
                end

                promotor_oi = species_option_struct.promotor_oi;
                species_name_list = species_option_struct.species_name_list; 
                conc_idx = species_option_struct.conc_idx; 
    
                % Sanity check: simulated data contains all species 
                num_species = length(species_name_list);
                sample_simData = Data{1};
                if ~isequal(size(sample_simData,2),num_species)
                    error('Mismatch in number of species in model')
                end

            % Extract datasets to be plotted 
            prom_idx = find(strcmp(promotor_name_list,promotor_oi)); 
                % No empty
            simulated_time_no_empty = Time{(prom_idx - 1) * 4 * length(conc_vec) + conc_idx}; 
            simulated_data_no_empty = Data{(prom_idx - 1) * 4 * length(conc_vec) + conc_idx}; 
                % Empty 
            simulated_time_empty = Time{(prom_idx - 1) * 4 * length(conc_vec) + 1 * length(conc_vec) + conc_idx}; 
            simulated_data_empty = Data{(prom_idx - 1) * 4 * length(conc_vec) + 1 * length(conc_vec) + conc_idx}; 
                % Empty T7 
            simulated_time_emptyT7 = Time{(prom_idx - 1) * 4 * length(conc_vec) + 2 * length(conc_vec) + conc_idx}; 
            simulated_data_emptyT7 = Data{(prom_idx - 1) * 4 * length(conc_vec) + 2 * length(conc_vec) + conc_idx}; 
                % Empty sigma 70 
            simulated_time_empty_sigma70 = Time{(prom_idx - 1) * 4 * length(conc_vec) + 3 * length(conc_vec) + conc_idx}; 
            simulated_data_empty_sigma70 = Data{(prom_idx - 1) * 4 * length(conc_vec) + 3 * length(conc_vec) + conc_idx}; 

            % Reorganize data 
                % Define species & resource of interest
            species_oi_list = {'sfGFP','broc','kanR','empty'};
            resource_oi_list = {'RNAP','t7RNAP','Ribo','RNase','AGTP','CUTP','AA','AGMP','CUMP','toxin'};
            bound_species_oi_name_list = {'sfGFP RNAP bound','broc RNAP bound','kanR RNAP bound','empty RNAP bound',...
                'sfGFP RNase bound','broc RNase bound','kanR RNase bound','empty RNase bound',...
                'sfGFP Ribo bound','broc Ribo bound','kanR Ribo bound','empty Ribo bound'};
            bound_species_name_list = {'RNAP','RNase','Ribo'}; 
            species_oi_total_conc_name_list = {'total sfGFP RNA','total broc RNA','total kanR RNA','total empty RNA'}; 
            
                % Preassign index list space for each named list (One plot
                % for each category) 
            resource_idx_list = cell(length(resource_oi_list),1); 
            bound_species_idx_list = cell(length(bound_species_oi_name_list),1);
            total_conc_idx_list = cell(length(species_oi_total_conc_name_list),1); 

            % Resource 
            for resource_species_idx = 1:length(resource_oi_list)
                resource_species_name = resource_oi_list{resource_species_idx};
                resource_idx_list{resource_species_idx} = find(strcmp(species_name_list,resource_species_name));
            end

            % Bound Species 

            for bound_species_idx = 1:length(bound_species_name_list)
                bound_species = bound_species_name_list{bound_species_idx}; 
                for species_oi_idx = 1:length(species_oi_list)
                    species_oi = species_oi_list{species_oi_idx}; 
                    bound_species_oi_idx = (bound_species_idx- 1) * length(species_oi_list) + species_oi_idx; 
                    selected_bound_species_idx_list = find(contains(species_name_list,species_oi) & contains(species_name_list,bound_species)); 
                    bound_species_idx_list{bound_species_oi_idx} = selected_bound_species_idx_list; 
                end
            end

            % Total Concentration 
            for total_conc_idx = 1:length(species_oi_total_conc_name_list)
                species_oi_name = species_oi_list{total_conc_idx};
                selected_species_oi_idx = contains(species_name_list,species_oi_name) & contains(species_name_list,'RNA') & ~contains(species_name_list,'RNase');
                total_conc_idx_list{total_conc_idx} = selected_species_oi_idx; 
            end

            % Plot 
                % Resource
            for resource_species_idx = 1:length(resource_oi_list)
                subplot(ceil(sqrt(length(resource_oi_list))),ceil(sqrt(length(resource_oi_list))),resource_species_idx)

                idx_oi_list = resource_idx_list{resource_species_idx}; 
                % Plot out 4 cases 
                plot(simulated_time_no_empty,sum(simulated_data_no_empty(:,idx_oi_list),2),'LineWidth',1.5,'Color',all_colors{1})
                hold on
                plot(simulated_time_empty,sum(simulated_data_empty(:,idx_oi_list),2),'LineWidth',1.5,'Color',all_colors{2})
                plot(simulated_time_emptyT7,sum(simulated_data_emptyT7(:,idx_oi_list),2),'LineWidth',1.5,'LineStyle','--','Color',all_colors{4})
                plot(simulated_time_empty_sigma70,sum(simulated_data_empty_sigma70(:,idx_oi_list),2),'LineWidth',1.5,'LineStyle','-.','Color',all_colors{3})

                species_name = resource_oi_list{resource_species_idx}; 
                species_name_for_plot = strrep(species_name,'_',' ');
                title(species_name_for_plot)
            end
            legend('No Empty','Empty','Empty T7','Empty sigma70')
            sgtitle('Resource')

                % Bound species
            figure; 
            for bound_species_idx = 1:length(bound_species_oi_name_list)
                subplot(ceil(sqrt(length(bound_species_oi_name_list))),ceil(sqrt(length(bound_species_oi_name_list))),bound_species_idx)

                idx_oi_list = bound_species_idx_list{bound_species_idx}; 

                % Plot out 4 cases 
                if ~isempty(idx_oi_list) 
                    data_flag = (max(sum(simulated_data_no_empty(:,idx_oi_list),2)) - min(sum(simulated_data_no_empty(:,idx_oi_list),2)) > 1e-08 | ...
                        max(sum(simulated_data_empty(:,idx_oi_list),2)) - min(sum(simulated_data_empty(:,idx_oi_list),2)) > 1e-08 | ...
                        max(sum(simulated_data_emptyT7(:,idx_oi_list),2)) - min(sum(simulated_data_emptyT7(:,idx_oi_list),2)) > 1e-08 | ...
                        max(sum(simulated_data_empty_sigma70(:,idx_oi_list),2)) - min(sum(simulated_data_empty_sigma70(:,idx_oi_list),2)) > 1e-08 )& ...
                        (max(abs(sum(simulated_data_no_empty(:,idx_oi_list),2))) > 1e-05 | max(abs(sum(simulated_data_empty(:,idx_oi_list),2))) > 1e-05 | ...
                        max(abs(sum(simulated_data_emptyT7(:,idx_oi_list),2))) > 1e-05 | max(abs(sum(simulated_data_empty_sigma70(:,idx_oi_list),2))) > 1e-05); 
                    if data_flag % Only plot if the value is not constant across time and has values 
                        plot(simulated_time_no_empty,sum(simulated_data_no_empty(:,idx_oi_list),2),'LineWidth',1.5,'Color',all_colors{1})
                        hold on
                        plot(simulated_time_empty,sum(simulated_data_empty(:,idx_oi_list),2),'LineWidth',1.5,'Color',all_colors{2})
                        plot(simulated_time_emptyT7,sum(simulated_data_emptyT7(:,idx_oi_list),2),'LineWidth',1.5,'LineStyle','--','Color',all_colors{4})
                        plot(simulated_time_empty_sigma70,sum(simulated_data_empty_sigma70(:,idx_oi_list),2),'LineWidth',1.5,'LineStyle','-.','Color',all_colors{3})
                    end
                end
                species_name = bound_species_oi_name_list{bound_species_idx}; 
                species_name_for_plot = strrep(species_name,'_',' ');
                title(species_name_for_plot)
            end
            legend('No Empty','Empty','Empty T7','Empty sigma70')
            sgtitle('Bound Species')

                % Total concentration 
            figure; 
            for total_conc_idx = 1:length(species_oi_total_conc_name_list)
                subplot(ceil(sqrt(length(species_oi_total_conc_name_list))),ceil(sqrt(length(species_oi_total_conc_name_list))),total_conc_idx)

                idx_oi_list = total_conc_idx_list{total_conc_idx}; 

                % Plot out 4 cases 
                if ~isempty(idx_oi_list) 
                    data_flag = (max(sum(simulated_data_no_empty(:,idx_oi_list),2)) - min(sum(simulated_data_no_empty(:,idx_oi_list),2)) > 1e-08 | ...
                        max(sum(simulated_data_empty(:,idx_oi_list),2)) - min(sum(simulated_data_empty(:,idx_oi_list),2)) > 1e-08 | ...
                        max(sum(simulated_data_emptyT7(:,idx_oi_list),2)) - min(sum(simulated_data_emptyT7(:,idx_oi_list),2)) > 1e-08 | ...
                        max(sum(simulated_data_empty_sigma70(:,idx_oi_list),2)) - min(sum(simulated_data_empty_sigma70(:,idx_oi_list),2)) > 1e-08) & ...
                        (max(abs(sum(simulated_data_no_empty(:,idx_oi_list),2))) > 1e-05 | max(abs(sum(simulated_data_empty(:,idx_oi_list),2))) > 1e-05 | ...
                        max(abs(sum(simulated_data_emptyT7(:,idx_oi_list),2))) > 1e-05 | max(abs(sum(simulated_data_empty_sigma70(:,idx_oi_list),2))) > 1e-05); 
                    if data_flag % Only plot if the value is not constant across time and has values 
                        plot(simulated_time_no_empty,sum(simulated_data_no_empty(:,idx_oi_list),2),'LineWidth',1.5,'Color',all_colors{1})
                        hold on
                        plot(simulated_time_empty,sum(simulated_data_empty(:,idx_oi_list),2),'LineWidth',1.5,'Color',all_colors{2})
                        plot(simulated_time_emptyT7,sum(simulated_data_emptyT7(:,idx_oi_list),2),'LineWidth',1.5,'LineStyle','--','Color',all_colors{4})
                        plot(simulated_time_empty_sigma70,sum(simulated_data_empty_sigma70(:,idx_oi_list),2),'LineWidth',1.5,'LineStyle','-.','Color',all_colors{3})
                    end
                end
                species_name = species_oi_total_conc_name_list{total_conc_idx}; 
                species_name_for_plot = strrep(species_name,'_',' ');
                title(species_name_for_plot)
            end
            legend('No Empty','Empty','Empty T7','Empty sigma70')
            sgtitle('Total Concentrations')

    end


end