timetotal = 43;
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

for tcount = 1:timetotal
   % clf
   if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    MCDS = read_MultiCellDS_xml(K);
    
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    deltax = abs(MCDS.mesh.X(1,1,k)-MCDS.mesh.X(1,2,k));
    deltay = abs(MCDS.mesh.Y(1,1,k)-MCDS.mesh.Y(2,1,k));
    
    virion(tcount) = sum(sum(MCDS.continuum_variables(1).data(:,:,k)))*deltax*deltay*20;%virion
    IFN(tcount) = sum(sum(MCDS.continuum_variables(2).data(:,:,k)))*deltax*deltay*20;%debris
    cytokine(tcount) = sum(sum(MCDS.continuum_variables(3).data(:,:,k)))*deltax*deltay*20;%proinflamcytokine virion
    chemokine(tcount) = sum(sum(MCDS.continuum_variables(4).data(:,:,k)))*deltax*deltay*20;%chemokine
    debris(tcount) = sum(sum(MCDS.continuum_variables(5).data(:,:,k)))*deltax*deltay*20;%debris
    ROS(tcount) = sum(sum(MCDS.continuum_variables(6).data(:,:,k)))*deltax*deltay*20;%debris
    
    total_cells = length(MCDS.discrete_cells.live_cells);
    
    CD8(tcount) = length(find( MCDS.discrete_cells.metadata.type == 3)); %CD8    
    neutrophils(tcount) = length(find( MCDS.discrete_cells.metadata.type == 5)); %neutrophils
    DC(tcount) = length(find( MCDS.discrete_cells.metadata.type == 6)); %DCs
    CD4(tcount) = length(find( MCDS.discrete_cells.metadata.type == 7)); %CD4
    
    inactivated_immune = find(MCDS.discrete_cells.custom.activated_immune_cell==0);
    activated_immune = find(MCDS.discrete_cells.custom.activated_immune_cell==1);
    
    macrophagesinactive(tcount) = length(intersect(find( MCDS.discrete_cells.metadata.type == 4),inactivated_immune)); %macs inactive
    macrophagesactive(tcount) = length(intersect(find( MCDS.discrete_cells.metadata.type == 4),activated_immune)); %macs active

    dead(tcount) = length(intersect(MCDS.discrete_cells.dead_cells,find(MCDS.discrete_cells.metadata.type == 1)));
    infected(tcount) = length(unique([find(MCDS.discrete_cells.custom.Vnuc>=1),MCDS.discrete_cells.dead_cells]));
    uninfected(tcount) = length(find( MCDS.discrete_cells.metadata.type == 1))-dead(tcount)-infected(tcount);
          
 end
 
 save('simulated_data.mat','virion','IFN','cytokine','chemokine','debris','ROS','CD8','neutrophils','DC','CD4','macrophagesinactive','macrophagesactive','dead','infected','uninfected')
