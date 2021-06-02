function solt = model_simulator_RK(parameters,initial_conditions)

solt = ode15s(@modemode,[0 12],initial_conditions); 

function dydt = modemode(t,y,Z)	
t
% Variables
VEx     = y(1);
VAtthi  = y(2);
VAttlo  = y(3);
VEn     = y(4);
Vcytp 	= y(5);
Vnucp   = y(6);
RC      = y(7);
RV      = y(8);	
RCRdRp  = y(9);
RVRdRp  = y(10);
Cp      = y(11);
VnucPM1 = y(12);
VcytPM1 = y(13);
RM1     = y(14);
RM2     = y(15);
RM3     = y(16);
RM4     = y(17);
RM5     = y(18);
RM6     = y(19);
RM7     = y(20);
RM8     = y(21);
PPB1    = y(22);
PPB2    = y(23);
PPA     = y(24);
PRdRp   = y(25);
PNP     = y(26);
PM1     = y(27);
PNEP    = y(28);
PHA     = y(29);
PNA     = y(30);
PM2     = y(31);
VRel    = y(32);
    
% Parameters
DRib = parameters.DRib;
		
FSp17 = parameters.FSp17;
FSp18 = parameters.FSp18;
	
kEn = parameters.kEn;
kEqhi = parameters.kEqhi;
kEqlo = parameters.kEqlo;
kFus = parameters.kFus;
kImp = parameters.kImp;
kExp = parameters.kExp;
kRdRp = parameters.kRdRp;
kRel = parameters.kRel;
	
kBindNP = parameters.kBindNP;
kBindM1 = parameters.kBindM1;
kBindRdRp = parameters.kBindRdRp;

kDegRnp = parameters.kDegRnp;
kDegR = parameters.kDegR;
kDegRRdRp = parameters.kDegRRdRp;
kDegM = parameters.kDegM;
	
kSynC = parameters.kSynC;
kSynV = parameters.kSynV;	
kSynM = parameters.kSynM;	
kSynP = parameters.kSynP;	
		
L1 = parameters.L1;
L2 = parameters.L2;
L3 = parameters.L3;
L4 = parameters.L4;
L5 = parameters.L5;
L6 = parameters.L6;
L7 = parameters.L7;
L8 = parameters.L8;
	
LV = parameters.LV;
		
KVRel = parameters.KVRel;
		
NNucNP = parameters.NNucNP;
NNucM1 = parameters.NNucM1;
NNucNEP = parameters.NNucNEP;
NPM1 = parameters.NPM1;
NPM2 = parameters.NPM2;
NPRdRp = parameters.NPRdRp;
NPHA = parameters.NPHA;
NPNP = parameters.NPNP;
NPNA = parameters.NPNA;
NPNEP = parameters.NPNEP;


kAtthi = parameters.kAtthi;
kAttlo = parameters.kAttlo;
FFus = parameters.FFus;
Btothi = parameters.Btothi;
Btotlo = parameters.Btotlo;

kDishi = kAtthi/kEqhi;
kDislo = kAttlo/kEqlo;

kDegVen = (1-FFus)/FFus*kFus;


%% calculate variables

Bhi = Btothi - VAtthi;
Blo = Btotlo - VAttlo;

% binding model
dVEx = kDishi*VAtthi+kDislo*VAttlo-(kAtthi*Bhi+kAttlo*Blo)*VEx;
dVAtthi = kAtthi*Bhi*VEx-(kDishi+kEn)*VAtthi;
dVAttlo = kAttlo*Blo*VEx-(kDislo+kEn)*VAttlo;
dVEn = kEn*(VAtthi+VAttlo)-(kFus+kDegVen)*VEn;

% replication model
rRel = kRel*VcytPM1*(PRdRp/(PRdRp+KVRel*NPRdRp))*(PHA/(PHA+KVRel*NPHA))*(PNP/(PNP+KVRel*NPNP))*(PNA/(PNA+KVRel*NPNA))*(PM1/(PM1+KVRel*NPM1))*(PM2/(PM2+KVRel*NPM2))*(PNEP/(PNEP+KVRel*NPNEP));
	
dVcytp   = 8*kFus*VEn-kImp*Vcytp;
dVnucp   = kImp*Vcytp+kBindNP*PNP*RVRdRp-(kBindM1*PM1+kDegRnp)*Vnucp;
dRC      = kSynC*Vnucp-kBindRdRp*PRdRp*RC-kDegR*RC;
dRV      = kSynV*Cp-kBindRdRp*PRdRp*RV-kDegR*RV;	
dRCRdRp  = kBindRdRp*PRdRp*RC-kBindNP*PNP*RCRdRp-kDegRRdRp*RCRdRp;
dRVRdRp  = kBindRdRp*PRdRp*RV-kBindNP*PNP*RVRdRp-kDegRRdRp*RVRdRp;
dCp      = kBindNP*PNP*RCRdRp-kDegRnp*Cp;
dVnucPM1 = kBindM1*PM1*Vnucp-(kExp*PNEP+kDegRnp)*VnucPM1;
dVcytPM1 = kExp*PNEP*VnucPM1-8*rRel-kDegRnp*VcytPM1;
dRM1     = (kSynM/L1)*(Vnucp/8)-kDegM*RM1;
dRM2     = (kSynM/L2)*(Vnucp/8)-kDegM*RM2;
dRM3     = (kSynM/L3)*(Vnucp/8)-kDegM*RM3;
dRM4     = (kSynM/L4)*(Vnucp/8)-kDegM*RM4;
dRM5     = (kSynM/L5)*(Vnucp/8)-kDegM*RM5;
dRM6     = (kSynM/L6)*(Vnucp/8)-kDegM*RM6;
dRM7     = (kSynM/L7)*(Vnucp/8)-kDegM*RM7;
dRM8     = (kSynM/L8)*(Vnucp/8)-kDegM*RM8;
dPPB1    = kSynP/DRib*RM2-kRdRp*PPB1*PPB2*PPA;
dPPB2    = kSynP/DRib*RM1-kRdRp*PPB1*PPB2*PPA;
dPPA     = kSynP/DRib*RM3-kRdRp*PPB1*PPB2*PPA;
dPRdRp   = kRdRp*PPB1*PPB2*PPA-kBindRdRp*PRdRp*(RV+RC)-(NPRdRp-8)*rRel;
dPNP     = kSynP/DRib*RM5-LV/NNucNP*kBindNP*PNP*(RVRdRp+RCRdRp)-(NPNP-8*LV/NNucNP)*rRel;
dPM1     = kSynP/DRib*(1-FSp17)*RM7-LV/NNucM1*kBindM1*PM1*Vnucp-(NPM1-8*LV/NNucM1)*rRel;
dPNEP    = kSynP/DRib*FSp18*RM8-LV/NNucNEP*kExp*PNEP*VnucPM1-(NPNEP-8*LV/NNucNEP)*rRel;
dPHA     = kSynP/DRib*RM4-NPHA*rRel;
dPNA     = kSynP/DRib*RM6-NPNA*rRel;
dPM2     = kSynP/DRib*FSp17*RM7-NPM2*rRel;
dVRel    = rRel;
   
dydt= [dVEx;dVAtthi;dVAttlo;dVEn;dVcytp;dVnucp;dRC;dRV;dRCRdRp;dRVRdRp;dCp;dVnucPM1;dVcytPM1;...
    dRM1;dRM2;dRM3;dRM4;dRM5;dRM6;dRM7;dRM8;dPPB1;dPPB2;dPPA;dPRdRp;dPNP;dPM1;dPNEP;dPHA;dPNA;dPM2;dVRel];
end
end