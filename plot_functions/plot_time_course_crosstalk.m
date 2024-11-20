function plot_time_course_crosstalk(sim_Time,sim_Data,exp_Time,exp_Data)

    % Plot time-course data in the crosstalk format

    promotor_name_list = {'T7_strong','T7_weak','sigma70_strong','sigma70_weak'};
    conc_vec = [0.5,1,2.5,5,10,15,30];
    all_colors = {'g',[0.5,0.5,0.5],'r','b'};
    num_promotors = length(sim_Data)/28; 
    figure;
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
                simulated_time = sim_Time{data_idx}; 
                simulated_data = sim_Data{data_idx}; 
        
                plot(simulated_time,simulated_data,'LineWidth',1.5,'LineStyle','--','Color',all_colors{empty_idx})
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
    line5 = plot(nan,nan,'Color','r','LineStyle','-','DisplayName','empty T7 experimental');
    line6 = plot(nan,nan,'Color','r','LineStyle','--','DisplayName','empty T7 simulated');
    line7 = plot(nan,nan,'Color','b','LineStyle','-','DisplayName','empty sigma70 experimental');
    line8 = plot(nan,nan,'Color','b','LineStyle','--','DisplayName','empty sigma70 simulated');
    if exist('exp_Time','var') && exist('exp_Data','var')
        legend([line1,line2,line3,line4,line5,line6,line7,line8])
    else
        legend([line2,line4,line6,line8])
    end




end