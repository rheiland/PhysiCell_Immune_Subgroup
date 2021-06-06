#include "./receptor_dynamics.h" 

using namespace PhysiCell; 

std::string receptor_model_version = "0.4.0"; 

Submodel_Information receptor_dynamics_info; 

void simple_receptor_dynamics_model_setup(void)
{

	// set version 
	receptor_dynamics_info.name = "receptor dynamics"; 
	receptor_dynamics_info.version = receptor_model_version; 
	// set functions 
	receptor_dynamics_info.main_function = simple_receptor_dynamics_main_model; 
	receptor_dynamics_info.phenotype_function = NULL; // pushed into the "main" model  
	receptor_dynamics_info.mechanics_function = NULL; 	
	
	// what microenvironment variables do you need 
	receptor_dynamics_info.microenvironment_variables.push_back( "virion" ); 	
	receptor_dynamics_info.microenvironment_variables.push_back( "interferon 1" ); 	
	receptor_dynamics_info.microenvironment_variables.push_back( "VTEST" ); 		
	
	// what custom data do I need? 
	receptor_dynamics_info.cell_variables.push_back( "VAtthi" ); 
	receptor_dynamics_info.cell_variables.push_back( "VAttlo" ); 
	receptor_dynamics_info.cell_variables.push_back( "Bhi" ); 	
	receptor_dynamics_info.cell_variables.push_back( "Blo" ); 	
	receptor_dynamics_info.cell_variables.push_back( "VEn" ); 	
	
	// submodel_registry.register_model( receptor_dynamics_info ); 
	receptor_dynamics_info.register_model(); 
	
	return;
}


void simple_receptor_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	if( phenotype.death.dead == true )
	{ return; } 


	static int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
		// if not lung epithelium, do nothing 
	if( pCell->type != lung_epithelial_type )
	{ return; } 
	
	//****************************************************************
	static int vtest_external = microenvironment.find_density_index( "VTEST" ); 
	
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
	double rho = pCell->nearest_density_vector()[vtest_external];
	double m = 1;//pCell->phenotype.molecular.internalized_total_substrates[vtest_external];
	double mhalf = 10;
	double rhomax = 100;
	
	if(rho>1/8000)
	{
		if(rho<rhomax/Vvoxel)
		{
			pCell->phenotype.secretion.uptake_rates[vtest_external] = parameters.doubles("uEvirus")*(mhalf/(m/Vvoxel+mhalf/Vvoxel));
		}
		else
		{
			pCell->phenotype.secretion.uptake_rates[vtest_external] = parameters.doubles("uEvirus")*(rhomax/Vvoxel/rho)*(mhalf/(m/Vvoxel+mhalf/Vvoxel));
		}
	}
	else 
	{pCell->phenotype.secretion.uptake_rates[vtest_external]=0;}
		
	

	/*
			
		// Determine the binding/unbinding rates from the model of Heldt and Laske to determine cell uptake/secretion rates
		
		// drawing variables for uptake from Heldt and Laske models 
		//double VEx = pCell->nearest_density_vector()[nV_external]*Vvoxel;//pCell->custom_data["VEx"];
		double VEx = pCell->nearest_density_vector()[vtest_external]*Vvoxel;//pCell->custom_data["VEx"];
		
		double VAtthi = pCell->custom_data["VAtthi"];
		double VAttlo = pCell->custom_data["VAttlo"];
		double Bhi = pCell->custom_data["Bhi"];
		double Blo = pCell->custom_data["Blo"];
		double VEn = pCell->custom_data["VEn"];
		double Vnuc = pCell->custom_data["Vnuc"];
		double VRel = pCell->custom_data["VRel"];
		
		// drawing parameters for uptake from Heldt and Laske models
		static double kAtthi = parameters.doubles("kAtthi");
		static double kAttlo = parameters.doubles("kAttlo");
		static double kDishi = parameters.doubles("kAtthi")/parameters.doubles("kEqhi");
		static double kDislo = parameters.doubles("kAttlo")/parameters.doubles("kEqlo");
		static double kEn = parameters.doubles("kEn");
		static double Btothi = parameters.doubles("Btothi");
		static double Btotlo = parameters.doubles("Btotlo");
		
		static double kRel = parameters.doubles("kRel");
		static double FFus = parameters.doubles("FFus");
		static double kFus = parameters.doubles("kFus");
		static double kDegVen = (1-FFus)/FFus*kFus;

		// evaluating Heldt and Laske's ODE models
		
		if(pCell->custom_data["antiviral_state"]>0.5) // cell is in an antiviral state =  no more uptake
		{
			
			// evaluating Heldt and Laske's ODE models
			pCell->custom_data["VAtthi"] += (-kDishi*VAtthi-kEn*VAtthi)*dt;
			pCell->custom_data["VAttlo"] += (-kDislo*VAttlo-kEn*VAttlo)*dt;
			pCell->custom_data["Bhi"] += Btothi-VAtthi;
			pCell->custom_data["Blo"] += Btotlo-VAttlo;
			pCell->custom_data["VEn"] += (kEn*(VAtthi+VAttlo)-(kFus-kDegVen)*VEn)*dt;												
			pCell->custom_data["VRel"] += (kRel*Vnuc)*dt;
			
			pCell->phenotype.secretion.secretion_rates[vtest_external] = (kDislo*VAttlo+kDishi*VAtthi+kRel*Vnuc)/Vvoxel; // does this need a dt or virions in the first two terms?
			pCell->phenotype.secretion.uptake_rates[vtest_external] = 0;
			
		}
		else if(Vnuc>1) // cell already replicating virus = no more uptake
		{
			// evaluating Heldt and Laske's ODE models
			pCell->custom_data["VAtthi"] += (-kDishi*VAtthi-kEn*VAtthi)*dt;
			pCell->custom_data["VAttlo"] += (-kDislo*VAttlo-kEn*VAttlo)*dt;
			pCell->custom_data["Bhi"] += Btothi-VAtthi;
			pCell->custom_data["Blo"] += Btotlo-VAttlo;
			pCell->custom_data["VEn"] += (kEn*(VAtthi+VAttlo)-(kFus-kDegVen)*VEn)*dt;										
			pCell->custom_data["VRel"] += (kRel*Vnuc)*dt;
			
			pCell->phenotype.secretion.secretion_rates[vtest_external] = (kDislo*VAttlo+kDishi*VAtthi+kRel*Vnuc)/Vvoxel; // does this need a dt or virions in the first two terms?
			pCell->phenotype.secretion.uptake_rates[vtest_external] = 0;
			
		}
		else if(VEx<0.9) // not enough virions outside to uptake
		{
			
			VEx = 0;
			pCell->custom_data["VAtthi"] += (-kDishi*VAtthi-kEn*VAtthi)*dt;
			pCell->custom_data["VAttlo"] += (-kDislo*VAttlo-kEn*VAttlo)*dt;
			pCell->custom_data["Bhi"] += Btothi-VAtthi;
			pCell->custom_data["Blo"] += Btotlo-VAttlo;
			pCell->custom_data["VEn"] += (kEn*(VAtthi+VAttlo)-(kFus-kDegVen)*VEn)*dt;
			pCell->custom_data["VRel"] += (kRel*Vnuc)*dt;
			
			pCell->phenotype.secretion.secretion_rates[vtest_external] = (kDislo*VAttlo+kDishi*VAtthi+kRel*Vnuc)/Vvoxel; // does this need a dt or virions in the first two terms?
			pCell->phenotype.secretion.uptake_rates[vtest_external] = 0;
			
			
		}
		else
		{		
			pCell->custom_data["VAtthi"] += (kAtthi*Bhi*VEx-kDishi*VAtthi-kEn*VAtthi)*dt;
			pCell->custom_data["VAttlo"] += (kAttlo*Blo*VEx-kDislo*VAttlo-kEn*VAttlo)*dt;
			pCell->custom_data["Bhi"] += Btothi-VAtthi;
			pCell->custom_data["Blo"] += Btotlo-VAttlo;
			pCell->custom_data["VEn"] += (kEn*(VAtthi+VAttlo)-(kFus-kDegVen)*VEn)*dt;
			pCell->custom_data["VRel"] += (kRel*Vnuc)*dt;
			
			pCell->phenotype.secretion.secretion_rates[vtest_external] = (kDislo*VAttlo+kDishi*VAtthi+kRel*Vnuc)/Vvoxel; // does this need a dt or virions in the first two terms?
			pCell->phenotype.secretion.uptake_rates[vtest_external] = (kAtthi*Bhi+kAttlo*Blo);
				
		}*/
		
		
	return;
}

void simple_receptor_dynamics_main_model( double dt )
{
	
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false )
		{ 
			simple_receptor_dynamics_model( pC, pC->phenotype , dt ); 
		}
	}
	
	return; 
}
