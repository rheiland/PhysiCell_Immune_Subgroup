#include "./internal_viral_dynamics.h" 

using namespace PhysiCell; 

std::string internal_virus_replication_version = "0.4.0"; 

Submodel_Information internal_viral_dynamics_info; 

void simple_internal_virus_model_setup(void)
{
	// set version
	internal_viral_dynamics_info.name = "internal viral replication dynamics"; 
	internal_viral_dynamics_info.version = internal_virus_replication_version; 
		// set functions 
	internal_viral_dynamics_info.main_function = NULL; 
	internal_viral_dynamics_info.phenotype_function = simple_internal_virus_model; 
	internal_viral_dynamics_info.mechanics_function = NULL; 
		// what microenvironment variables do I need? 

		// what custom data do I need? 
	internal_viral_dynamics_info.microenvironment_variables.push_back( "assembled virion" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "VEn" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "Virions_internalized" ); 	
	

	internal_viral_dynamics_info.cell_variables.push_back( "Vcytp" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "Vnucp" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RC" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RV" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RCRdRp" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RVRdRp" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "Cp" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "VnucPM1" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "VcytPM1" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "RM1" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RM2" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RM3" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RM4" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RM5" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RM6" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RM7" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "RM8" ); 		
	internal_viral_dynamics_info.cell_variables.push_back( "PPB1" ); 		
	internal_viral_dynamics_info.cell_variables.push_back( "PPB2" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "PPA" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "PRdRp" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "PNP" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "PM1" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "PNEP" );
	internal_viral_dynamics_info.cell_variables.push_back( "PHA" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "PNA" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "PM2" ); 
	internal_viral_dynamics_info.cell_variables.push_back( "VRel" ); 	
	internal_viral_dynamics_info.cell_variables.push_back( "rRel" ); 			 			

	// submodel_registry.register_model( internal_viral_dynamics_info ); 
	internal_viral_dynamics_info.register_model();
	return;
}


void simple_internal_virus_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{ return; } 

	static int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
		// if not lung epithelium, do nothing 
	if( pCell->type != lung_epithelial_type )
	{ return; } 

	// Simple virus replication model
	static int IFN_index = microenvironment.find_density_index( "interferon 1" );
	static int nV_external = microenvironment.find_density_index( "virion" );
	double gamma = parameters.doubles("gamma");
	double IFN_internal = pCell->nearest_density_vector()[IFN_index];
	double IC_50_IFN = parameters.doubles("IC_50_IFN");
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
	
	double m_i = pCell->phenotype.molecular.internalized_total_substrates[nV_external];

	if(IFN_internal<0)
	{IFN_internal = 0;}
	double gamma_IFN = gamma;//*(1-IFN_internal/(IC_50_IFN+IFN_internal));
		
	intracellular_replication_model(pCell, phenotype, dt);
	
	//std::cout<<pCell->custom_data["VEn"]<<std::endl;	
	// Changing viral replication model to match heldt uptake model
	//std::cout<<pCell->custom_data["VEn"]<<std::endl;
	//if(pCell->custom_data["VEn"]>1e-5)
	//{
		//std::cout<<pCell->custom_data["VEn"]<<std::endl;
	//	pCell->custom_data["VEn"] += gamma_IFN*dt;
	//}
	//pCell->custom_data["Virions_internalized"] = pCell->custom_data["VEn"];
	
	return;
}

void intracellular_replication_model(  Cell* pCell, Phenotype& phenotype, double dt )
{
		// drawing variables for uptake from Heldt and Laske models 
	double VEn = pCell->custom_data["VEn"];
	double Vcytp = pCell->custom_data["Vcytp"];
	double Vnucp = pCell->custom_data["Vnucp"];
	double RC = pCell->custom_data["RC"];
	double RV = pCell->custom_data["RV"];	
	double RCRdRp = pCell->custom_data["RCRdRp"];
	double RVRdRp = pCell->custom_data["RVRdRp"];
	double Cp = pCell->custom_data["Cp"];
	double VnucPM1 = pCell->custom_data["VnucPM1"];
	double VcytPM1 = pCell->custom_data["VcytPM1"];
	double RM1 = pCell->custom_data["RM1"];
	double RM2 = pCell->custom_data["RM2"];
	double RM3 = pCell->custom_data["RM3"];
	double RM4 = pCell->custom_data["RM4"];
	double RM5 = pCell->custom_data["RM5"];
	double RM6 = pCell->custom_data["RM6"];
	double RM7 = pCell->custom_data["RM7"];
	double RM8 = pCell->custom_data["RM8"];
	double PPB1 = pCell->custom_data["PPB1"];
	double PPB2 = pCell->custom_data["PPB2"];
	double PPA = pCell->custom_data["PPA"];
	double PRdRp = pCell->custom_data["PRdRp"];
	double PNP = pCell->custom_data["PNP"];
	double PM1 = pCell->custom_data["PM1"];
	double PNEP = pCell->custom_data["PNEP"];
	double PHA = pCell->custom_data["PHA"];
	double PNA = pCell->custom_data["PNA"];
	double PM2 = pCell->custom_data["PM2"];
	double VRel = pCell->custom_data["VRel"];
		
	// drawing parameters for uptake from Heldt and Laske models
	
	static double DRib = parameters.doubles("DRib");
		
	static double FSp17 = parameters.doubles("FSp17");
	static double FSp18 = parameters.doubles("FSp18");
	
	static double kEn = parameters.doubles("kEn");
	static double kEqhi = parameters.doubles("kEqhi");
	static double kEqlo = parameters.doubles("kEqlo");
	static double kFus = parameters.doubles("kFus");
	static double kImp = parameters.doubles("kImp");
	static double kExp = parameters.doubles("kExp");
	static double kRdRp = parameters.doubles("kRdRp");
	static double kRel = parameters.doubles("kRel");
	
	static double kBindNP = parameters.doubles("kBindNP");
	static double kBindM1 = parameters.doubles("kBindM1");
	static double kBindRdRp = parameters.doubles("kBindRdRp");
	
	static double kDegRnp = parameters.doubles("kDegRnp");
	static double kDegR = parameters.doubles("kDegR");
	static double kDegRRdRp = parameters.doubles("kDegRRdRp");
	static double kDegM = parameters.doubles("kDegM");
	
	static double kSynC = parameters.doubles("kSynC");
	static double kSynV = parameters.doubles("kSynV");	
	static double kSynM = parameters.doubles("kSynM");	
	static double kSynP = parameters.doubles("kSynP");	
		
	static double L1 = parameters.doubles("L1");
	static double L2 = parameters.doubles("L2");
	static double L3 = parameters.doubles("L3");
	static double L4 = parameters.doubles("L4");
	static double L5 = parameters.doubles("L5");
	static double L6 = parameters.doubles("L6");
	static double L7 = parameters.doubles("L7");
	static double L8 = parameters.doubles("L8");
	
	static double LV = parameters.doubles("LV");
		
	static double KVRel = parameters.doubles("KVRel");
		
	static double NNucNP = parameters.doubles("NNucNP");
	static double NNucM1 = parameters.doubles("NNucM1");
	static double NNucNEP = parameters.doubles("NNucNEP");
	static double NPM1 = parameters.doubles("NPM1");
	static double NPM2 = parameters.doubles("NPM2");
	static double NPRdRp = parameters.doubles("NPRdRp");
	static double NPHA = parameters.doubles("NPHA");
	static double NPNP = parameters.doubles("NPNP");
	static double NPNA = parameters.doubles("NPNA");
	static double NPNEP = parameters.doubles("NPNEP");
	
		
	// Equations
	
	double rRel = kRel*VcytPM1;
	rRel *=(PRdRp/(PRdRp+KVRel*NPRdRp));
	rRel *=(PHA/(PHA+KVRel*NPHA));
	rRel *=(PNP/(PNP+KVRel*NPNP));
	rRel *=(PNA/(PNA+KVRel*NPNA));
	rRel *=(PM1/(PM1+KVRel*NPM1));
	rRel *=(PM2/(PM2+KVRel*NPM2));
	rRel *=(PNEP/(PNEP+KVRel*NPNEP));
	
	pCell->custom_data["rRel"] = rRel;
	
	pCell->custom_data["Vcytp"] += (8*kFus*VEn-kImp*Vcytp)*dt;
	pCell->custom_data["Vnucp"] += (kImp*Vcytp+kBindNP*PNP*RVRdRp-(kBindM1*PM1+kDegRnp)*Vnucp)*dt;
	pCell->custom_data["RC"] += (kSynC*Vnucp-kBindRdRp*PRdRp*RC-kDegR*RC)*dt;
	pCell->custom_data["RV"] += (kSynV*Cp-kBindRdRp*PRdRp*RV-kDegR*RV)*dt;	
	pCell->custom_data["RCRdRp"] += (kBindRdRp*PRdRp*RC-kBindNP*PNP*RCRdRp-kDegRRdRp*RCRdRp)*dt;
	pCell->custom_data["RVRdRp"] += (kBindRdRp*PRdRp*RV-kBindNP*PNP*RVRdRp-kDegRRdRp*RVRdRp)*dt;
	pCell->custom_data["Cp"] += (kBindNP*PNP*RCRdRp-kDegRnp*Cp)*dt;
	pCell->custom_data["VnucPM1"] += (kBindM1*PM1*Vnucp-(kExp*PNEP+kDegRnp)*VnucPM1)*dt;
	pCell->custom_data["VcytPM1"] += (kExp*PNEP*VnucPM1-8*rRel-kDegRnp*VcytPM1)*dt;
	pCell->custom_data["RM1"] += ((kSynM/L1)*(Vnucp/8)-kDegM*RM1)*dt;
	pCell->custom_data["RM2"] += ((kSynM/L2)*(Vnucp/8)-kDegM*RM2)*dt;
	pCell->custom_data["RM3"] += ((kSynM/L3)*(Vnucp/8)-kDegM*RM3)*dt;
	pCell->custom_data["RM4"] += ((kSynM/L4)*(Vnucp/8)-kDegM*RM4)*dt;
	pCell->custom_data["RM5"] += ((kSynM/L5)*(Vnucp/8)-kDegM*RM5)*dt;
	pCell->custom_data["RM6"] += ((kSynM/L6)*(Vnucp/8)-kDegM*RM6)*dt;
	pCell->custom_data["RM7"] += ((kSynM/L7)*(Vnucp/8)-kDegM*RM7)*dt;
	pCell->custom_data["RM8"] += ((kSynM/L8)*(Vnucp/8)-kDegM*RM8)*dt;
	pCell->custom_data["PPB1"] += (kSynP/DRib*RM2-kRdRp*PPB1*PPB2*PPA)*dt;
	pCell->custom_data["PPB2"] += (kSynP/DRib*RM1-kRdRp*PPB1*PPB2*PPA)*dt;
	pCell->custom_data["PPA"] += (kSynP/DRib*RM3-kRdRp*PPB1*PPB2*PPA)*dt;
	pCell->custom_data["PRdRp"] += (kRdRp*PPB1*PPB2*PPA-kBindRdRp*PRdRp*(RV+RC)-(NPRdRp-8)*rRel)*dt;
	pCell->custom_data["PNP"] += (kSynP/DRib*RM5-LV/NNucNP*kBindNP*PNP*(RVRdRp+RCRdRp)-(NPNP-8*LV/NNucNP)*rRel)*dt;
	pCell->custom_data["PM1"] += (kSynP/DRib*(1-FSp17)*RM7-LV/NNucM1*kBindM1*PM1*Vnucp-(NPM1-8*LV/NNucM1)*rRel)*dt;
	pCell->custom_data["PNEP"] += (kSynP/DRib*FSp18*RM8-LV/NNucNEP*kExp*PNEP*VnucPM1-(NPNEP-8*LV/NNucNEP)*rRel)*dt;
	pCell->custom_data["PHA"] += (kSynP/DRib*RM4-NPHA*rRel)*dt;
	pCell->custom_data["PNA"] += (kSynP/DRib*RM6-NPNA*rRel)*dt;
	pCell->custom_data["PM2"] += (kSynP/DRib*FSp17*RM7-NPM2*rRel)*dt;
	pCell->custom_data["VRel"] += (rRel)*dt;

	
	if(pCell->custom_data["VRel"]>1e-6)
	{std::cout<<"VRel: "<<pCell->custom_data["VRel"]<<std::endl;}

	return;
	
}
