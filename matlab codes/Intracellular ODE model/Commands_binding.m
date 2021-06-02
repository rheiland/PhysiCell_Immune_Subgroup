

% initialise variables
VEx     = 5*20*20*20;
VAtthi  = 0; 
VAttlo  = 0; 
VEn     = 0;
VV      = 0;
VRel     = 0;

%set parameters	
parameters.kEn = 4.8/60;
parameters.kEqhi = 1.13*1e-2;
parameters.kEqlo = 8.33*1e-5;
parameters.kFus = 3.21/60;
parameters.kAtthi = 8.09*1e-2/60;
parameters.kAttlo = 4.55*1e-4/60; %xxxxxxx 4.55*1e-4;%/60;
parameters.FFus = 0.51;
parameters.Btothi = 150;
parameters.Btotlo = 1000;
parameters.kRel = 3.7*1e-3/60;



parameters.eps_gam = 1e4;
parameters.gam = 1.39/60;


dt = 0.01;
iter = 72000;
timegrid = linspace(0,dt*iter,iter)/60;


% Recreating Figure 2

sol_mat = [VEx;VAtthi;VAttlo;VEn;VV;VRel];

for i = 1:iter
   
   val(:,i) = sol_mat;
   sol_mat_updated = model_simulator_binding(dt,sol_mat,parameters);	
   sol_mat = sol_mat_updated;
    
   
end

figure
hold on
plot(timegrid,val(1,:)')
plot(timegrid,val(2,:)')
plot(timegrid,val(3,:)')
plot(timegrid,val(4,:)')
yyaxis right
plot(timegrid,val(5,:)')
plot(timegrid,val(6,:)')
legend('VEx','VAtthi','VAttlo','VEn','VV','VRel')
title('Figure 2')

VVoxel = 20*20*20;

figure
hold on 
yyaxis left
plot(timegrid,val(1,:)'/VVoxel)
yyaxis right 
plot(timegrid,val(6,:)'/VVoxel)


%%

dt = 0.01;
iter = 72000

for i = 1:iter
   
   val(:,i) = sol_mat;
   sol_mat_updated = model_simulator_discritised_2(dt,sol_mat,parameters);	
   sol_mat = sol_mat_updated;
    
   
end

timegrid = linspace(0,dt*iter,iter);%/60;

figure
hold on 
plot(timegrid,val')

