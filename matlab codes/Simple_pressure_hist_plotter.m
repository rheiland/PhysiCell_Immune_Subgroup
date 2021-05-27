
MCDS = read_MultiCellDS_xml('output00000002.xml');
figure
hist(MCDS.discrete_cells.custom.displacement_stor,100)
xlabel('Simple pressure')
ylabel('Frequency')
set(gca,'FontSize',18)
grid on

[h bins] = hist(MCDS.discrete_cells.custom.displacement_stor,100);

[freq_max, max_loc] = max(h);
av_simp_pressure = bins(max_loc)

bounadry_cells_xU = find(MCDS.discrete_cells.state.position(:,1)>374);
bounadry_cells_xL = find(MCDS.discrete_cells.state.position(:,1)<-374);
bounadry_cells_yU = find(MCDS.discrete_cells.state.position(:,2)>364);
bounadry_cells_yL = find(MCDS.discrete_cells.state.position(:,2)<-374);

figure
hold on 
plot(MCDS.discrete_cells.state.position(bounadry_cells_xU,1),MCDS.discrete_cells.state.position(bounadry_cells_xU,2),'o')
plot(MCDS.discrete_cells.state.position(bounadry_cells_xL,1),MCDS.discrete_cells.state.position(bounadry_cells_xL,2),'o')
plot(MCDS.discrete_cells.state.position(bounadry_cells_yU,1),MCDS.discrete_cells.state.position(bounadry_cells_yU,2),'o')
plot(MCDS.discrete_cells.state.position(bounadry_cells_yL,1),MCDS.discrete_cells.state.position(bounadry_cells_yL,2),'o')

cells_on_boundary  = [bounadry_cells_xU;bounadry_cells_xL;bounadry_cells_yU;bounadry_cells_yL];

figure
hist(MCDS.discrete_cells.custom.displacement_stor(cells_on_boundary),100)
xlabel('Simple pressure')
ylabel('Frequency')
set(gca,'FontSize',18)
grid on

Simple_pressure = MCDS.discrete_cells.custom.displacement_stor;
Simple_pressure(cells_on_boundary) = [];

figure
hist(Simple_pressure,100)
xlabel('Simple pressure')
ylabel('Frequency')
set(gca,'FontSize',18)
grid on
