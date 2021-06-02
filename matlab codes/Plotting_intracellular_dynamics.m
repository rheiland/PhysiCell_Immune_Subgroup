timetotal = 48;
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
  
    VEx(tcount) = MCDS.discrete_cells.custom.VEx(3242);
    VAtthi(tcount) = MCDS.discrete_cells.custom.VAtthi(3242);
    VAttlo(tcount) = MCDS.discrete_cells.custom.VAttlo(3242);
    VEn(tcount) = MCDS.discrete_cells.custom.VEn(3242);
    Vnuc(tcount) = MCDS.discrete_cells.custom.Vnuc(3242);
    VRel(tcount) = MCDS.discrete_cells.custom.VRel(3242);

    
 end

 figure
 hold on 
 plot(VEx)
 title('VEx')
 
 figure
 hold on 
 plot(VAtthi)
 title('VAtthi')
  
 figure
 hold on 
 plot(VAttlo)
 title('VAttlo')

 figure
 hold on 
 plot(VEn)
 title('VEn')
 
 figure
 hold on 
 plot(Vnuc)
 title('Vnuc')
 
 figure
 hold on 
 plot(VRel)
 title('VRel')
 