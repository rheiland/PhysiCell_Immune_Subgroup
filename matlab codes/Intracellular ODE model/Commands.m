

% initialise variables
VEx     = 10;
VAtthi  = 0; 
VAttlo  = 0; 
VEn     = 0;
Vcytp 	= 0;
Vnucp   = 0;
RC      = 0;
RV      = 0;	
RCRdRp  = 0;
RVRdRp  = 0;
Cp      = 0;
VnucPM1 = 0;
VcytPM1 = 0;
RM1     = 0;
RM2     = 0;
RM3     = 0;
RM4     = 0;
RM5     = 0;
RM6     = 0;
RM7     = 0;
RM8     = 0;
PPB1    = 0;
PPB2    = 0;
PPA     = 0;
PRdRp   = 0;
PNP     = 0;
PM1     = 0;
PNEP    = 0;
PHA     = 0;
PNA     = 0;
PM2     = 0;
VRel    = 0;

rRel = 0;


%set parameters
parameters.DRib = 160;
		
parameters.FSp17 = 0.02;
parameters.FSp18 = 0.125;
	
parameters.kEn = 4.8;%/60;
parameters.kEqhi = 1.13*1e-2;
parameters.kEqlo = 8.33*1e-5;
parameters.kFus = 3.21;%/60;
parameters.kImp = 6;%/60;
parameters.kExp = 1e-6;%/60;
parameters.kRdRp = 1;%/60;
parameters.kRel = 3.7*1e-3;%/60;
	
parameters.kBindNP = 3.01*1e-4;%/60;
parameters.kBindM1 = 1.39*1e-6;%/60;
parameters.kBindRdRp = 1;%/60;

parameters.kDegRnp = 0.09;%/60;
parameters.kDegR = 36.36;%/60;
parameters.kDegRRdRp = 4.25;%/60;
parameters.kDegM = 0.33;%/60;
	
parameters.kSynC = 1.38;%/60;
parameters.kSynV = 13.86;%/60;	
parameters.kSynM = 2.5*1e5;%/60;	xxxxx 2.5*1e5;%/60;
parameters.kSynP = 64800;%/60;	
		
parameters.L1 = 2320;
parameters.L2 = 2320;
parameters.L3 = 2211;
parameters.L4 = 1757;
parameters.L5 = 1540;
parameters.L6 = 1392;
parameters.L7 = 1005;
parameters.L8 = 868;
	
parameters.LV = 1700;
		
parameters.KVRel = 10;
		
parameters.NNucNP = 24;
parameters.NNucM1 = 200;
parameters.NNucNEP = 1700;
parameters.NPM1 = 3000;
parameters.NPM2 = 40;
parameters.NPRdRp = 45;
parameters.NPHA = 500;
parameters.NPNP = 1000;
parameters.NPNA = 100;
parameters.NPNEP = 165;

parameters.kAtthi = 8.09*1e-2;%/60;
parameters.kAttlo = 4.55*1e-4;%/60; xxxxxxx 4.55*1e-4;%/60;
parameters.FFus = 0.51;
parameters.Btothi = 150;
parameters.Btotlo = 1000;

dt = 1.6667e-04;%0.01;
iter = 72000;
timegrid = linspace(0,dt*iter,iter);%/60;


% Recreating Figure 2

VEx     = 10;
VAtthi  = 0; 
VAttlo  = 0; 

sol_mat = [VEx;VAtthi;VAttlo;VEn;Vcytp;Vnucp;RC;RV;RCRdRp;RVRdRp;Cp;VnucPM1;VcytPM1;...
    RM1;RM2;RM3;RM4;RM5;RM6;RM7;RM8;PPB1;PPB2;PPA;PRdRp;PNP;PM1;PNEP;PHA;PNA;PM2;VRel];

for i = 1:iter
   
   val(:,i) = sol_mat;
   sol_mat_updated = model_simulator_discritised(dt,sol_mat,parameters);	
   sol_mat = sol_mat_updated;
    
   
end

figure
yyaxis left
hold on
plot(timegrid,val(1,:)')
plot(timegrid,val(4,:)')
yyaxis right
plot(timegrid,val(6,:)')
xlim([0 2])
legend('VEx','VEn','Vnucp')
title('Figure 2')


%% compared with ODE45 sim

% initialise variables
VEx     = 0;
VAtthi  = 10; 
VAttlo  = 0; 
VEn     = 0;
Vcytp 	= 0;
Vnucp   = 0;
RC      = 0;
RV      = 0;	
RCRdRp  = 0;
RVRdRp  = 0;
Cp      = 0;
VnucPM1 = 0;
VcytPM1 = 0;
RM1     = 0;
RM2     = 0;
RM3     = 0;
RM4     = 0;
RM5     = 0;
RM6     = 0;
RM7     = 0;
RM8     = 0;
PPB1    = 0;
PPB2    = 0;
PPA     = 0;
PRdRp   = 0;
PNP     = 0;
PM1     = 0;
PNEP    = 0;
PHA     = 0;
PNA     = 0;
PM2     = 0;
VRel    = 0;

rRel = 0;


%set parameters
parameters.DRib = 160;
		
parameters.FSp17 = 0.02;
parameters.FSp18 = 0.125;
	
parameters.kEn = 4.8;%/60;
parameters.kEqhi = 1.13*1e-2;
parameters.kEqlo = 8.33*1e-5;
parameters.kFus = 3.21;%/60;
parameters.kImp = 6;%/60;
parameters.kExp = 1e-6;%/60;
parameters.kRdRp = 1;%/60;
parameters.kRel = 3.7*1e-3;%/60;
	
parameters.kBindNP = 3.01*1e-4;%/60;
parameters.kBindM1 = 1.39*1e-6;%/60;
parameters.kBindRdRp = 1;%/60;

parameters.kDegRnp = 0.09;%/60;
parameters.kDegR = 36.36;%/60;
parameters.kDegRRdRp = 4.25;%/60;
parameters.kDegM = 0.33;%/60;
	
parameters.kSynC = 1.38;%/60;
parameters.kSynV = 13.86;%/60;	
parameters.kSynM = 2.5*1e5;%/60;	%XXXXX
parameters.kSynP = 64800;%/60;	
		
parameters.L1 = 2320;
parameters.L2 = 2320;
parameters.L3 = 2211;
parameters.L4 = 1757;
parameters.L5 = 1540;
parameters.L6 = 1392;
parameters.L7 = 1005;
parameters.L8 = 868;
	
parameters.LV = 1700;
		
parameters.KVRel = 10;
		
parameters.NNucNP = 24;
parameters.NNucM1 = 200;
parameters.NNucNEP = 1700;
parameters.NPM1 = 3000;
parameters.NPM2 = 40;
parameters.NPRdRp = 45;
parameters.NPHA = 500;
parameters.NPNP = 1000;
parameters.NPNA = 100;
parameters.NPNEP = 165;

parameters.kAtthi = 8.09*1e-2;%/60;
parameters.kAttlo = 4.55*1e-4;%/60;  % XXXXXX
parameters.FFus = 0.51;
parameters.Btothi = 150;
parameters.Btotlo = 1000;


parameters.kEqhi = 1e8;
VEx =100;
VAtthi = 0;

sol_mat = [VEx;VAtthi;VAttlo;VEn;Vcytp;Vnucp;RC;RV;RCRdRp;RVRdRp;Cp;VnucPM1;VcytPM1;...
    RM1;RM2;RM3;RM4;RM5;RM6;RM7;RM8;PPB1;PPB2;PPA;PRdRp;PNP;PM1;PNEP;PHA;PNA;PM2;VRel];

initial_conditions = [VEx;VAtthi;VAttlo;VEn;Vcytp;Vnucp;RC;RV;RCRdRp;RVRdRp;Cp;VnucPM1;VcytPM1;...
    RM1;RM2;RM3;RM4;RM5;RM6;RM7;RM8;PPB1;PPB2;PPA;PRdRp;PNP;PM1;PNEP;PHA;PNA;PM2;VRel];

sol = model_simulator_RK(parameters,initial_conditions);

figure
yyaxis left
hold on
plot(sol.x,sol.y(1,:)','LineWidth',2)
plot(sol.x,sol.y(4,:)','LineWidth',2)
yyaxis right
plot(sol.x,sol.y(6,:)','LineWidth',2)
%xlim([0 3])
legend('VEx','VEn','Vcytp','Vnucp')
title('Figure 2')


figure
yyaxis left
plot(sol.x,sol.y(8,:)','LineWidth',2)
yyaxis right
plot(sol.x,sol.y(10,:)','LineWidth',2)
legend('RV','RVRdRp')

figure
yyaxis left
plot(sol.x,sol.y(6,:)','LineWidth',2)
yyaxis right
hold on
plot(sol.x,sol.y(27,:)','LineWidth',2)
plot(sol.x,sol.y(26,:)','LineWidth',2)
plot(sol.x,sol.y(25,:)','LineWidth',2)
legend('Vcytp','PM1','PNP','PRdRp')


figure
yyaxis left
plot(sol.x,sol.y(13,:)','LineWidth',2)
yyaxis right
hold on
plot(sol.x,sol.y(32,:)','LineWidth',2)
legend('VpcytM1','VRel')

%%

% initialise variables
VEx     = 0;
VAtthi  = 10; 
VAttlo  = 0; 
VEn     = 0;
Vcytp 	= 0;
Vnucp   = 0;
RC      = 0;
RV      = 0;	
RCRdRp  = 0;
RVRdRp  = 0;
Cp      = 0;
VnucPM1 = 0;
VcytPM1 = 0;
RM1     = 0;
RM2     = 0;
RM3     = 0;
RM4     = 0;
RM5     = 0;
RM6     = 0;
RM7     = 0;
RM8     = 0;
PPB1    = 0;
PPB2    = 0;
PPA     = 0;
PRdRp   = 0;
PNP     = 0;
PM1     = 0;
PNEP    = 0;
PHA     = 0;
PNA     = 0;
PM2     = 0;
VRel    = 0;

rRel = 0;


%set parameters
parameters.DRib = 160;
		
parameters.FSp17 = 0.02;
parameters.FSp18 = 0.125;
	
parameters.kEn = 4.8/60;
parameters.kEqhi = 1.13*1e-2;
parameters.kEqlo = 8.33*1e-5;
parameters.kFus = 3.21/60;
parameters.kImp = 6/60;
parameters.kExp = 1e-6/60;
parameters.kRdRp = 1/60;
parameters.kRel = 3.7*1e-3/60;
	
parameters.kBindNP = 3.01*1e-4/60;
parameters.kBindM1 = 1.39*1e-6/60;
parameters.kBindRdRp = 1/60;

parameters.kDegRnp = 0.09/60;
parameters.kDegR = 36.36/60;
parameters.kDegRRdRp = 4.25/60;
parameters.kDegM = 0.33/60;
	
parameters.kSynC = 1.38/60;
parameters.kSynV = 13.86/60;	
parameters.kSynM = 2.5*1e5/60;	%XXXXX
parameters.kSynP = 64800/60;	
		
parameters.L1 = 2320;
parameters.L2 = 2320;
parameters.L3 = 2211;
parameters.L4 = 1757;
parameters.L5 = 1540;
parameters.L6 = 1392;
parameters.L7 = 1005;
parameters.L8 = 868;
	
parameters.LV = 1700;
		
parameters.KVRel = 10;
		
parameters.NNucNP = 24;
parameters.NNucM1 = 200;
parameters.NNucNEP = 1700;
parameters.NPM1 = 3000;
parameters.NPM2 = 40;
parameters.NPRdRp = 45;
parameters.NPHA = 500;
parameters.NPNP = 1000;
parameters.NPNA = 100;
parameters.NPNEP = 165;

parameters.kAtthi = 8.09*1e-2/60;
parameters.kAttlo = 4.55*1e-4/60;  % XXXXXX
parameters.FFus = 0.51;
parameters.Btothi = 150;
parameters.Btotlo = 1000;


parameters.kEqhi = 1e8;
VEx =100;
VAtthi = 0;

%-------
sol_mat = [VEx;VAtthi;VAttlo;VEn;Vcytp;Vnucp;RC;RV;RCRdRp;RVRdRp;Cp;VnucPM1;VcytPM1;...
    RM1;RM2;RM3;RM4;RM5;RM6;RM7;RM8;PPB1;PPB2;PPA;PRdRp;PNP;PM1;PNEP;PHA;PNA;PM2;VRel];

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

