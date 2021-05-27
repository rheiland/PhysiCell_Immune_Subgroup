timetotal = 240;
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

 for tcount = 1:timetotal+1
   % clf
   if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    MCDS = read_MultiCellDS_xml(K);
    
  ind1(tcount) = length(find( MCDS.discrete_cells.metadata.type == 3)); %CD8
  inactivated_immune = find(MCDS.discrete_cells.custom.activated_immune_cell==0);
  activated_immune = find(MCDS.discrete_cells.custom.activated_immune_cell==1);
  ind2in(tcount) = length(intersect(find( MCDS.discrete_cells.metadata.type == 4),inactivated_immune)); %macs inactive
  ind2ac(tcount) = length(intersect(find( MCDS.discrete_cells.metadata.type == 4),activated_immune)); %macs active
  ind3(tcount) = length(find( MCDS.discrete_cells.metadata.type == 5)); %neutrophils
  ind4(tcount) = length(find( MCDS.discrete_cells.metadata.type == 6)); %DCs
  ind5(tcount) = length(find( MCDS.discrete_cells.metadata.type == 7)); %CD4
  
    dead_cells(tcount) = length(MCDS.discrete_cells.dead_cells);
    infected_cells(tcount) = length(unique([find(MCDS.discrete_cells.custom.Virions_internalized>=1),MCDS.discrete_cells.dead_cells]));
    uninfected = find(MCDS.discrete_cells.custom.Virions_internalized<1);
    live_target_cells(tcount) = length(find( MCDS.discrete_cells.metadata.type == 1))-dead_cells(tcount)-infected_cells(tcount);
   
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    deltax = abs(MCDS.mesh.X(1,1,k)-MCDS.mesh.X(1,2,k));
    deltay = abs(MCDS.mesh.Y(1,1,k)-MCDS.mesh.Y(2,1,k));
    virion(tcount) = sum(sum(MCDS.continuum_variables(1).data(:,:,k)))*deltax*deltay*20;%virion
   assembled_virion(tcount) = sum(sum(MCDS.continuum_variables(2).data(:,:,k)))*deltax*deltay*20;%assembled virion
   pro_inflam_cytokine(tcount) = sum(sum(MCDS.continuum_variables(4).data(:,:,k)))*deltax*deltay*20;%proinflamcytokine virion
   chemokine(tcount) = sum(sum(MCDS.continuum_variables(5).data(:,:,k)))*deltax*deltay*20;%chemokine
   debris(tcount) = sum(sum(MCDS.continuum_variables(6).data(:,:,k)))*deltax*deltay*20;%debris

    IFN(tcount) = sum(sum(MCDS.continuum_variables(3).data(:,:,k)))*deltax*deltay*20;%debris
    ROS(tcount) = sum(sum(MCDS.continuum_variables(7).data(:,:,k)))*deltax*deltay*20;%debris

    
 end

 hours_to_days = linspace(0,timetotal/(24*1),timetotal+1);
 
 figure
 hold on 
 plot(hours_to_days,live_target_cells,'LineWidth',2)
 plot(hours_to_days,infected_cells,'--','LineWidth',2)
 plot(hours_to_days,dead_cells,':','LineWidth',2)
 ylabel('Number of cells')
% ylim([0 3000])
 yyaxis right
 plot(hours_to_days,virion','-.','LineWidth',2)
 ylabel('Total virions')
 legend('Uninfected cells','Infected cells','Dead cells','Virions')
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 saveas(gcf,'F1.png')
 

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

 figure
 hold on
 plot(hours_to_days,ind2in,'Color',[0.04 0.31 0.49],'LineWidth',2)
 plot(hours_to_days,ind2ac,'--','Color',[.3 .75 .93],'LineWidth',2)
 plot(hours_to_days,ind3,':','Color',[.47 .67 .19],'LineWidth',2)
 legend('Macrophages (inactive)','Macrophages (active)','Neutrophils','CD8 T cells')
 xlabel('Time (days)')
 ylabel('Number of cells')
 set(gca,'FontSize',15)
 saveas(gcf,'F2.png')
 
figure
 hold on
 plot(hours_to_days,ind4,'Color',[1 0 0],'LineWidth',2)
 plot(hours_to_days,ind1,'--','Color',[1 0.07 0.65],'LineWidth',2)
 plot(hours_to_days,ind5,':','Color',[.64 .08 .18],'LineWidth',2)
  legend('DCs','CD8 T cells','CD4 T cells')
 xlabel('Time (days)')
 ylabel('Number of cells')
 set(gca,'FontSize',15)
 saveas(gcf,'F3.png')
 
 figure
 hold on 
 plot(hours_to_days,pro_inflam_cytokine,'Color',[1 0.41 0.16],'LineWidth',2)
 plot(hours_to_days,chemokine,'--','Color',[0.12 0.64 0.54],'LineWidth',2)
 plot(hours_to_days,debris,':','Color',[.73 0.46 0.9],'LineWidth',2)
 plot(hours_to_days,IFN,':','LineWidth',2)
 plot(hours_to_days,ROS,':','LineWidth',2)
 ylabel('Total substrates')
 legend('Pro-inflammatory cytokine','Chemokine','Debris','IFN','ROS')
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 saveas(gcf,'F4.png')
 
 STOP 
 figure
 hold on 
 plot(hours_to_days,live_target_cells,'LineWidth',2)
 plot(hours_to_days,infected_cells,'--','LineWidth',2)
 plot(hours_to_days,dead_cells,':','LineWidth',2)
 ylim([0 3000])
 ylabel('Number of cells')
 yyaxis right
 plot(hours_to_days,virion,'Color',[0 0.45 0.74],'LineWidth',2)
 ylabel('Total virions')
 legend('Live cells','Infected cells','Dead cells','Virion')
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 
 %%
 
 
load('small_tissue_RNA_threshold_1')
 load('small_tissue_RNA_threshold_5')
 load('small_tissue_RNA_threshold_100')
 
 hours_to_days = linspace(0,timetotal/(24*2),timetotal);
 
 figure
 hold on 
 plot(hours_to_days,live_target_cells_thresh_1,'LineWidth',2)
 plot(hours_to_days,live_target_cells_thresh_5,'LineWidth',2)
 plot(hours_to_days,live_target_cells_thresh_100,'LineWidth',2)
 ylabel('Number of live cells')
 ylim([0 3000])
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 
 figure
 hold on 
 plot(hours_to_days,infected_cells_thresh_1,'LineWidth',2)
 plot(hours_to_days,infected_cells_thresh_5,'LineWidth',2)
 plot(hours_to_days,infected_cells_thresh_100,'LineWidth',2)
 ylabel('Number of infected cells')
 ylim([0 3000])
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 
  
 figure
 hold on 
 plot(hours_to_days,dead_cells_thresh_1,'LineWidth',2)
 plot(hours_to_days,dead_cells_thresh_5,'LineWidth',2)
 plot(hours_to_days,dead_cells_thresh_100,'LineWidth',2)
 ylabel('Number of dead cells')
 ylim([0 3000])
  xlabel('Time (days)')
 set(gca,'FontSize',15)

 
 figure
 hold on
 plot(hours_to_days,virion_thresh_1','LineWidth',2)
 plot(hours_to_days,virion_thresh_5','LineWidth',2)
 plot(hours_to_days,virion_thresh_100','LineWidth',2)
 ylabel('Total virions')
 xlabel('Time (days)')
 set(gca,'FontSize',15)

 
 figure
 hold on 
 plot(hours_to_days,pro_inflam_cytokine_thresh_1,'LineWidth',2)
 plot(hours_to_days,pro_inflam_cytokine_thresh_5,'LineWidth',2)
 plot(hours_to_days,pro_inflam_cytokine_thresh_100,'LineWidth',2)
 ylabel('Total cytokine')
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 legend('s_{threshold} = 1','s_{threshold}=5','s_{threshold} = 100')
 
 
 
 
 figure
 hold on 
 plot(hours_to_days,live_target_cells,'LineWidth',2)
 plot(hours_to_days,infected_cells,'--','LineWidth',2)
 plot(hours_to_days,dead_cells,':','LineWidth',2)
 ylim([0 3000])
 ylabel('Number of cells')
 yyaxis right
 plot(hours_to_days,virion,'Color',[0 0.45 0.74],'LineWidth',2)
 ylabel('Total virions')
 legend('Live cells','Infected cells','Dead cells','Virion')
 xlabel('Time (days)')
 set(gca,'FontSize',15)
 
 
 
