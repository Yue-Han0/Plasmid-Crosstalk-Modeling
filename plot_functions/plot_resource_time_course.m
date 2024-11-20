function plot_resource_time_course(Time,Data,resource_oi,resource_names,legend_names,figure_title,save_path)

    % Plot out resource utilization for a set of data (designed for plasmid
    % dosage) 

    % Inputs 
        % Time: a cell array of nT * 1 time vectors 
        % Data: a cell array of nT * #tracked resources data vectors 
        % resource_oi: a cell array of tracked resources in model
        % resouce_names: a cell array of conventional names for tracked
        % resouces 

    figure; 
    for resource_idx = 1:length(resource_names)
        subplot(ceil(sqrt(length(resource_names))),ceil(sqrt(length(resource_names))),resource_idx)
        for conc_idx = 1:length(legend_names)
            Time_single = Time{conc_idx,1}; 
            Data_single = Data{conc_idx,1}; 
            if length(resource_oi) > length(resource_names) && isequal(resource_idx,length(resource_names))
                plot(Time_single,sum(Data_single(:,resource_idx:end),2),'LineWidth',1.5)
            else
                plot(Time_single,Data_single(:,resource_idx),'LineWidth',1.5)
            end
            hold on 
        end
        title(resource_names{resource_idx})

    end

    legend(legend_names)


    if exist("figure_title",'var')
        sgtitle(figure_title)
    end

    if exist('save_path','var')
        saveas(gcf,save_path,'png')
    end




end