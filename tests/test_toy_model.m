clear
clc
% Modified basic gene expression with T7 promotor expressing sfGFP 
% plasmid_conc_list = [0.1,1,5,10,15,30,75,100];
% plasmid_conc_list = [0.5,1]; 
% plasmid_conc_list = [0.5,1,2.5,5,10,15,30];
plasmid_conc_list= [0.5,5,30]; 
% plasmid_conc_list = 10 .* plasmid_conc_list; 
% empty_plasmid_conc = 15; 
empty_plasmid_conc = 15; 

% plasmid_conc_list = [50,100,200,500,1000,1500,3000];

% 
% t7RNAP_init_cond = 100; 
RNAP_init_cond = 750; 
RNase_init_cond = 10;
% Ribo_init_cond = 100; 
% % RNase_init_cond = 30000; 
% 
RNase_reporter_R = 3690; 
RNase_kanR_R = 3.69;
RNase_empty_R = 3.69; 
% % 
RNA_deg = 0.087664 / 10; 


% empty_plasmid_conc = 300; 
% set path 
currentpath = pwd;
addpath(genpath(sprintf('%s/txtlsim_buildacell',currentpath))); 
addpath(genpath(sprintf('%s/plasmid_crosstalk_config_files',currentpath))); 
rmpath(sprintf('%s/txtlsim_buildacell/components',currentpath));
% plasmid_conc_list = [1e-05,1e-04,1e-03,1e-02,0.1]; 
figure; 
for plasmid_conc_idx = 1:length(plasmid_conc_list)
    % Set up the standard TXTL tubes
    % These load up the RNAP, Ribosome and degradation enzyme concentrations
    tube1 = txtl_extract('plasmid_crosstalk');
    tube2 = txtl_buffer('plasmid_crosstalk');
    
    %% No empty 
    % Now set up a tube that will contain our DNA
    tube3 = txtl_newtube('no_empty_gene_expression');
    
    % Define the tube with only reporter plasmid 
    % dna_deGFP = txtl_add_dna(tube3, 'pJ23119(35)', 'utrGFP(57)', 'sfGFP(723)',plasmid_conc_list(plasmid_conc_idx), 'plasmid');	
    % dna_broc = txtl_add_dna(tube3, 'pT7(23)','utrbroc(57)', 'no_protein',plasmid_conc_list(plasmid_conc_idx), 'plasmid');	
    dna_broc = txtl_add_dna(tube3, 'pJ23119(35)', 'utrbroc(57)', 'no_protein',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    % dna_sfGFP = txtl_add_dna(tube3, 'pT7(23)', 'utrGFP(57)', 'sfGFP(723)',plasmid_conc_list(plasmid_conc_idx), 'plasmid');	
    dna_kan = txtl_add_dna(tube3,'pkanR(20)','utrkanR(30)','kanR(783)',plasmid_conc_list(plasmid_conc_idx),'plasmid');
    
    % Mix the contents of the individual tubes
    Mobj = txtl_combine([tube1, tube2, tube3]);
    [~] = txtl_runsim(Mobj,14*60*60);
    
    %% Empty 
    % Define another tube with reporter plasmid + kan 
    tube4 = txtl_newtube('empty_gene_expression');
    % dna_deGFP = txtl_add_dna(tube4, 'pJ23119(35)', 'utrGFP(57)', 'sfGFP(723)',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    % dna_broc = txtl_add_dna(tube4, 'pT7(23)', 'utrbroc(57)', 'no_protein',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    dna_broc = txtl_add_dna(tube4, 'pJ23119(35)', 'utrbroc(57)', 'no_protein',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    % dna_deGFP_2 = txtl_add_dna(tube4, 'pT7(23)', 'utrGFP(57)', 'sfGFP(723)',plasmid_conc_list(plasmid_conc_idx), 'plasmid');	
    dna_kan_2 = txtl_add_dna(tube4,'pkanR(20)','utrkanR(30)','kanR(783)',empty_plasmid_conc+plasmid_conc_list(plasmid_conc_idx),'plasmid');
    Mobj_empty = txtl_combine([tube1, tube2, tube4]); 
    [~] = txtl_runsim(Mobj_empty,14*60*60); 

    % 
    %% Empty T7 
    tube5 = txtl_newtube('emptyT7_gene_expression');
    
    % dna_deGFP = txtl_add_dna(tube5, 'pJ23119(35)', 'utrGFP(57)', 'sfGFP(723)',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    % dna_deGFP = txtl_add_dna(tube5, 'pT7(23)', 'utrGFP(57)', 'sfGFP(723)',plasmid_conc_list(plasmid_conc_idx), 'plasmid');	
    % dna_broc = txtl_add_dna(tube5, 'pT7(23)', 'utrbroc(57)', 'no_protein',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    dna_broc = txtl_add_dna(tube5, 'pJ23119(35)', 'utrbroc(57)', 'no_protein',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    dna_empty = txtl_add_dna(tube5, 'pT7(23)', 'utrempty(51)', 'no_protein',empty_plasmid_conc, 'plasmid');	
    dna_kan = txtl_add_dna(tube5,'pkanR(20)','utrkanR(30)','kanR(783)',empty_plasmid_conc + plasmid_conc_list(plasmid_conc_idx),'plasmid');
    Mobj_emptyT7 = txtl_combine([tube1, tube2, tube5]);
    [~] = txtl_runsim(Mobj_emptyT7,14*60*60); 

    %% Empty sigma70 
    tube6 = txtl_newtube('empty_sigma70_gene_expression');
    
    % dna_deGFP = txtl_add_dna(tube6, 'pJ23119(35)', 'utrGFP(57)', 'sfGFP(723)',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    % dna_broc = txtl_add_dna(tube6, 'pT7(23)', 'utrbroc(57)', 'no_protein',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    dna_broc = txtl_add_dna(tube6, 'pJ23119(35)', 'utrbroc(57)', 'no_protein',plasmid_conc_list(plasmid_conc_idx), 'plasmid');
    % dna_deGFP = txtl_add_dna(tube6, 'pT7(23)', 'utrGFP(57)', 'sfGFP(723)',plasmid_conc_list(plasmid_conc_idx), 'plasmid');	
    dna_empty = txtl_add_dna(tube6, 'pJ23119(35)', 'utrempty(51)', 'no_protein',empty_plasmid_conc, 'plasmid');	
    dna_kan = txtl_add_dna(tube6,'pkanR(20)','utrkanR(30)','kanR(783)',empty_plasmid_conc + plasmid_conc_list(plasmid_conc_idx),'plasmid');
    Mobj_empty_sigma70 = txtl_combine([tube1, tube2, tube6]);
    [~] = txtl_runsim(Mobj_empty_sigma70,14*60*60);

    % Need to set T7RNAP concentration to a nonzero value 
        % Modify initial amount field for Mobj.Species with the name
    % % 't7RNAP' 
    % set(sbioselect(Mobj.Species,'Name','t7RNAP'),'InitialAmount',t7RNAP_init_cond);
    % set(sbioselect(Mobj_empty.Species,'Name','t7RNAP'),'InitialAmount',t7RNAP_init_cond);
    % set(sbioselect(Mobj_emptyT7.Species,'Name','t7RNAP'),'InitialAmount',t7RNAP_init_cond);
    % set(sbioselect(Mobj_empty_sigma70.Species,'Name','t7RNAP'),'InitialAmount',t7RNAP_init_cond);
    % % RNAP
    % set(sbioselect(Mobj.Species,'Name','RNAP'),'InitialAmount',RNAP_init_cond);
    % set(sbioselect(Mobj_empty.Species,'Name','RNAP'),'InitialAmount',RNAP_init_cond);
    % set(sbioselect(Mobj_emptyT7.Species,'Name','RNAP'),'InitialAmount',RNAP_init_cond);
    % set(sbioselect(Mobj_empty_sigma70.Species,'Name','RNAP'),'InitialAmount',RNAP_init_cond);
    % % RNase
    % set(sbioselect(Mobj.Species,'Name','RNase'),'InitialAmount',RNase_init_cond);
    % set(sbioselect(Mobj_empty.Species,'Name','RNase'),'InitialAmount',RNase_init_cond);
    % set(sbioselect(Mobj_emptyT7.Species,'Name','RNase'),'InitialAmount',RNase_init_cond);
    % set(sbioselect(Mobj_empty_sigma70.Species,'Name','RNase'),'InitialAmount',RNase_init_cond);
    % % Ribosome 
    % set(sbioselect(Mobj.Species,'Name','Ribo'),'InitialAmount',Ribo_init_cond);
    % set(sbioselect(Mobj_empty.Species,'Name','Ribo'),'InitialAmount',Ribo_init_cond);
    % set(sbioselect(Mobj_emptyT7.Species,'Name','Ribo'),'InitialAmount',Ribo_init_cond);
    % set(sbioselect(Mobj_empty_sigma70.Species,'Name','Ribo'),'InitialAmount',Ribo_init_cond);

    % % RNase_binding_rxn_name = '[RNA utrbroc--no_protein] + RNase <-> [RNA utrbroc--no_protein:RNase]'; 
    % RNase_binding_rxn = Mobj.Reactions(25);
    % RNase_binding_param_R = sbioselect(RNase_binding_rxn.KineticLaw.Parameters,'Name','TXTL_RNAdeg_R');
    % set(RNase_binding_param_R,'Value',3.69e+03)
    % 
    % RNase_binding_rxn = Mobj_empty.Reactions(25);
    % RNase_binding_param_R = sbioselect(RNase_binding_rxn.KineticLaw.Parameters,'Name','TXTL_RNAdeg_R');
    % set(RNase_binding_param_R,'Value',3.69e+03)
    % 
    % RNase_binding_rxn = Mobj_emptyT7.Reactions(35);
    % RNase_binding_param_R = sbioselect(RNase_binding_rxn.KineticLaw.Parameters,'Name','TXTL_RNAdeg_R');
    % set(RNase_binding_param_R,'Value',3.69e+03)
    % % RNase_binding_rxn = Mobj_emptyT7.Reactions(19);
    % % RNase_binding_param_R = sbioselect(RNase_binding_rxn.KineticLaw.Parameters,'Name','TXTL_RNAdeg_R');
    % % set(RNase_binding_param_R,'Value',3.69e+02)
    % 
    % RNase_binding_rxn = Mobj_empty_sigma70.Reactions(35);
    % RNase_binding_param_R = sbioselect(RNase_binding_rxn.KineticLaw.Parameters,'Name','TXTL_RNAdeg_R');
    % set(RNase_binding_param_R,'Value',3.69e+03)
    % RNase_binding_rxn = Mobj_empty_sigma70.Reactions(19);
    % RNase_binding_param_R = sbioselect(RNase_binding_rxn.KineticLaw.Parameters,'Name','TXTL_RNAdeg_R');
    % set(RNase_binding_param_R,'Value',3.69e+02)


    % set(sbioselect(Mobj,'Type','Parameter','Name','AGTPdeg_time'),'Value',3000);
    % set(sbioselect(Mobj_empty,'Type','Parameter','Name','AGTPdeg_time'),'Value',3000);
    % set(sbioselect(Mobj_emptyT7,'Type','Parameter','Name','AGTPdeg_time'),'Value',3000);
    % set(sbioselect(Mobj_empty_sigma70,'Type','Parameter','Name','AGTPdeg_time'),'Value',3000);
    
    % Run a simulaton
    %
    % At this point, the entire experiment is set up and loaded into 'Mobj'.
    % So now we just use standard Simbiology and MATLAB commands to run
    % and plot our results
    %
    
    % Mobj = modify_RNase_binding_parameters(Mobj,RNase_reporter_R,RNase_kanR_R);
    % Mobj_empty = modify_RNase_binding_parameters(Mobj_empty,RNase_reporter_R,RNase_kanR_R);
    % Mobj_emptyT7 = modify_RNase_binding_parameters(Mobj_emptyT7,RNase_reporter_R,RNase_kanR_R,RNase_empty_R);
    % Mobj_empty_sigma70 = modify_RNase_binding_parameters(Mobj_empty_sigma70,RNase_reporter_R,RNase_kanR_R,RNase_empty_R);
    % 
    % Mobj = modify_RNA_deg_parameters(Mobj,RNA_deg);
    % Mobj_empty = modify_RNA_deg_parameters(Mobj_empty,RNA_deg);
    % Mobj_emptyT7 = modify_RNA_deg_parameters(Mobj_emptyT7,RNA_deg);
    % Mobj_empty_sigma70 = modify_RNA_deg_parameters(Mobj_empty_sigma70,RNA_deg);

    tic
    [simData] = txtl_runsim(Mobj,14*60*60);
    % temp_plot_resource_change(simData)
    toc
    
    tic 
    [simData_empty] = txtl_runsim(Mobj_empty,14*60*60); 
    % temp_plot_resource_change(simData_empty)
    toc

    tic
    [simData_emptyT7] = txtl_runsim(Mobj_emptyT7,14*60*60); 
    % temp_plot_resource_change(simData_emptyT7)
    % % txtl_plot(simData_emptyT7,Mobj_emptyT7)
    toc
    
    tic
    [simData_empty_sigma70] = txtl_runsim(Mobj_empty_sigma70,14*60*60); 
    % temp_plot_resource_change(simData_empty_sigma70)
    % txtl_plot(simData_empty_sigma70,Mobj_empty_sigma70)
    toc

    % figure(1);
    
    hold on 
    subplot(2,2,plasmid_conc_idx)
    % mRNA_idx = find(contains(simData.DataNames,'RNA utrGFP--sfGFP')); 
    % mRNA_idx = find(contains(simData.DataNames,'RNA utrGFP--sfGFP') & ~contains(simData.DataNames,'RNase')); 
    % protein_idx = find(strcmp(simData.DataNames,'protein sfGFP*'));
    mRNA_idx = find(strcmp(simData.DataNames,'RNA utrbroc--no_protein')); 
    % mRNA_idx = find(strcmp(simData.Data))
    % mRNA_idx = find(contains(simData.DataNames,'RNA utrGFP--sfGFP') & ~contains(simData.DataNames,'RNase')); 
    plot(simData.Time./3600,sum(simData.Data(:,mRNA_idx),2),'LineWidth',1.5,'Color','g')
    hold on 
    mRNA_idx = find(contains(simData_empty.DataNames,'RNA utrbroc--no_protein') & ~contains(simData_empty.DataNames,'RNase')); 
    plot(simData_empty.Time./3600,sum(simData_empty.Data(:,mRNA_idx),2),'LineWidth',1.5,'Color',[0.5,0.5,0.5],'Marker','+')
    % hold on 
    % mRNA_idx = find(contains(simData_emptyT7.DataNames,'RNA utrbroc--no_protein') & ~contains(simData_emptyT7.DataNames,'RNase')); 
    % plot(simData_emptyT7.Time,sum(simData_emptyT7.Data(:,mRNA_idx),2),'LineWidth',1.5,'Color','r','Marker','o')
    % hold on 
    % mRNA_idx = find(contains(simData_empty_sigma70.DataNames,'RNA utrbroc--no_protein') & ~contains(simData_empty_sigma70.DataNames,'RNase')); 
    % plot(simData_empty_sigma70.Time,sum(simData_empty_sigma70.Data(:,mRNA_idx),2),'LineWidth',1.5,'Color','b')
    % title(sprintf('%.4f nM reporter plasmid',plasmid_conc_list(plasmid_conc_idx)))
    hold off 
    % sgtitle('protein broc T7')
    sgtitle('mRNA broc sigma70')
    if isequal(plasmid_conc_idx,length(plasmid_conc_list))
        legend('Reporter only (no empty)','Reporter + kanR expression (empty)','empty T7','empty sigma70')
    end
    ylabel('Concentration(nM)','FontSize',14)
    xlabel('Time(hr)','FontSize',14)
    set(gca,'FontSize',14)
    legend('Baseline','w/ Empty Vector')

    % % 2024/03/13 edits: add resource plots 
    % figure; 
    % all_species_name_list = {Mobj.Species.Name}; 
    % mRNA_kanR_idx_list = contains(all_species_name_list,'RNA utrkanR--kanR') & ~contains(all_species_name_list,'RNase');
    % mRNA_kanR_names = all_species_name_list(mRNA_kanR_idx_list); 
    % resource_oi_name_list = [{'RNAP','RNase','AGTP','RNA utrbroc--no_protein'} mRNA_kanR_names];
    % resource_oi_formal_name_list = {'RNA Polymerase','Ribonuclease','AGTP','RNA reporter','RNA kanR'};
    % for species_idx = 1:5
    %     subplot(2,3,species_idx)
    % 
    %     % Get name for resource of interest 
    %     if ~isequal(species_idx,5)
    %         species_oi_name = resource_oi_name_list{species_idx};
    %     else
    %         species_oi_name = resource_oi_name_list(species_idx:end); 
    %     end
    % 
    %     % Get index for resource of interest 
    %     if iscell(species_oi_name)
    %         species_oi_idx = nan(1,length(species_oi_name)); 
    %         for i = 1:length(species_oi_name)
    %             species_oi_name_single = species_oi_name{i}; 
    %             species_oi_idx(i) = find(strcmp(all_species_name_list,species_oi_name_single)); 
    %         end
    %     else
    %         species_oi_idx = strcmp(all_species_name_list,species_oi_name); 
    %     end
    % 
    %     % Get data and plot
    %         % Baseline
    %     species_oi_data_baseline = sum(simData.Data(:,species_oi_idx),2); 
    %     timeVec_baseline = simData.Time; 
    %         % with empty vector 
    %     species_oi_data_empty = sum(simData_empty.Data(:,species_oi_idx),2); 
    %     timeVec_empty = simData_empty.Time; 
    %         % Plot
    %     plot(timeVec_baseline./3600,species_oi_data_baseline,'LineWidth',1.5,'Color','g')
    %     hold on 
    %     plot(timeVec_empty./3600,species_oi_data_empty,'LineWidth',1.5,'Color',[0.5,0.5,0.5])
    %     title(resource_oi_formal_name_list{species_idx})
    %     xlabel('Time (hr)','FontSize',14)
    %     ylabel('Concentration (nM)','FontSize',14)
    %     set(gca,'FontSize',14)
    % 
    % 
    % end
    

%     figure(2);
%     hold on 
%     subplot(3,3,plasmid_conc_idx)
%     protein_idx = find(strcmp(simData.DataNames,'protein sfGFP*')); 
%     plot(simData.time,simData.Data(:,protein_idx),'LineWidth',1.5,'Color','k')
%     hold on 
%     plot(simData_empty.time,simData_empty.Data(:,protein_idx),'LineWidth',1.5,'Color','r')
%     title(sprintf('%.1f nM reporter plasmid',plasmid_conc_list(plasmid_conc_idx)))
%     hold off
%     sgtitle('Protein sfGFP sigma70')
%     if isequal(plasmid_conc_idx,length(plasmid_conc_list))
%         legend('Reporter only (no empty)','Reporter + kanR expression (empty)')
%     end
    
%     txtl_plot(simData,Mobj);
%     txtl_plot(simData_empty,Mobj_empty);
% 
%     1 + 1 ==2 
end
% plot the result
% txtl_plot(simData,Mobj);

% Automatically use matlab mode in emacs (keep at end of file)
% Local variables:
% mode: matlab
% End:


