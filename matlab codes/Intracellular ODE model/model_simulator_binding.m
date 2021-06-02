function sol_mat_updated = model_simulator_discritised(dt,sol_mat,parameters)	

% Variables
VEx     = sol_mat(1);
VAtthi  = sol_mat(2);
VAttlo  = sol_mat(3);
VEn     = sol_mat(4);
VV      = sol_mat(5);
VRel    = sol_mat(6);

% Parameters
kEn = parameters.kEn;
kRel = parameters.kRel;
kEqhi = parameters.kEqhi;
kEqlo = parameters.kEqlo;
kFus = parameters.kFus;
kAtthi = parameters.kAtthi;
kAttlo = parameters.kAttlo;
FFus = parameters.FFus;
Btothi = parameters.Btothi;
Btotlo = parameters.Btotlo;
kDishi = kAtthi/kEqhi;
kDislo = kAttlo/kEqlo;

gam = parameters.gam;
eps_gam = parameters.eps_gam;

kDegVen = (1-FFus)/FFus*kFus;

%% calculate variables

Bhi = Btothi - VAtthi;
Blo = Btotlo - VAttlo;

% binding model
VEx = VEx + (kDishi*VAtthi+kDislo*VAttlo-(kAtthi*Bhi+kAttlo*Blo)*VEx)*dt;
VAtthi = VAtthi + (kAtthi*Bhi*VEx-(kDishi+kEn)*VAtthi)*dt;
VAttlo = VAttlo + (kAttlo*Blo*VEx-(kDislo+kEn)*VAttlo)*dt;
VEn = VEn + (kEn*(VAtthi+VAttlo)-(kFus+kDegVen)*VEn)*dt;
VV = VV+(kFus*VEn+gam*VV*(1-VV/(eps_gam))-kRel*VV)*dt;
VRel = VRel+(kRel*VV)*dt;
   
sol_mat_updated= [VEx;VAtthi;VAttlo;VEn;VV;VRel];
end