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
	
	// submodel_registry.register_model( receptor_dynamics_info ); 
	receptor_dynamics_info.register_model(); 
	
	return;
}

/*
void receptor_dynamics_model_setup( void )
{
		// set version 
	receptor_dynamics_info.name = "receptor dynamics"; 
	receptor_dynamics_info.version = receptor_model_version; 
		// set functions 
	receptor_dynamics_info.main_function = receptor_dynamics_main_model; 
	receptor_dynamics_info.phenotype_function = NULL; // pushed into the "main" model  
	receptor_dynamics_info.mechanics_function = NULL; 	
		// what microenvironment variables do you need 
	receptor_dynamics_info.microenvironment_variables.push_back( "virion" ); 		
	receptor_dynamics_info.microenvironment_variables.push_back( "assembled virion" ); 		
		// what cell variables and parameters do you need? 
	receptor_dynamics_info.cell_variables.push_back( "unbound_external_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "bound_external_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "unbound_internal_ACE2" ); 
	receptor_dynamics_info.cell_variables.push_back( "bound_internal_ACE2" ); 
	
	receptor_dynamics_info.cell_variables.push_back( "ACE2_binding_rate" ); 
	receptor_dynamics_info.cell_variables.push_back( "ACE2_endocytosis_rate" ); 
	receptor_dynamics_info.cell_variables.push_back( "ACE2_cargo_release_rate" ); 	
	receptor_dynamics_info.cell_variables.push_back( "ACE2_recycling_rate" ); 
	
	// submodel_registry.register_model( receptor_dynamics_info ); 
	receptor_dynamics_info.register_model(); 
	
	return; 
}*/

void simple_receptor_dynamics_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{ return; } 

	static int lung_epithelial_type = get_cell_definition( "lung epithelium" ).type; 
		// if not lung epithelium, do nothing 
	if( pCell->type != lung_epithelial_type )
	{ return; } 

	
	//****************************************************************
	
	//std::cout<<"here 1"<<std::endl;
	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int IFN_index = microenvironment.find_density_index("interferon 1");
	
	double rho_virus = pCell->nearest_density_vector()[nV_external];
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
	double m_i = pCell->phenotype.molecular.internalized_total_substrates[nV_external]*Vvoxel;//pCell->custom_data["intracellular_virus_amount"];//pCell->custom_data["virions"];
	double m_half = parameters.doubles("m_half");
	double rho_half = parameters.doubles("rho_half");
	
	double u_Evirus = parameters.doubles("u_Evirus");
	double IFN_internal = pCell->nearest_density_vector()[IFN_index];
	double IC_50_IFN = parameters.doubles("IC_50_IFN");
	double rho_max = parameters.doubles("rho_max");



	if(  IFN_internal<1e-16)
	{
		IFN_internal = 0;
	}
	double u_Evirus_IFN = u_Evirus*(1-IFN_internal/(IC_50_IFN+IFN_internal));
	
	
	// Determine the binding/unbinding rates from the model of Heldt and Laske to determine cell uptake/secretion rates
	
	// drawing variables for uptake from Heldt and Laske models 
	double VEx = pCell->custom_data["VEx"];
	double VAtthi = pCell->custom_data["VAtthi"];
	double VAttlo = pCell->custom_data["VAttlo"];
	double Bhi = pCell->custom_data["Bhi"];
	double Blo = pCell->custom_data["Blo"];
	double VEn = pCell->custom_data["VEn"];
	
	// drawing parameters for uptake from Heldt and Laske models
	double kAtthi = parameters.doubles("kAtthi");
	double kAttlo = parameters.doubles("kAttlo");
	double kDishi = parameters.doubles("kAtthi")/parameters.doubles("kEqhi");
	double kDislo = parameters.doubles("kAttlo")/parameters.doubles("kEqlo");
	double kEn = parameters.doubles("kEn");
	double Btothi = parameters.doubles("Btothi");
	double Btotlo = parameters.doubles("Btotlo");
	
	double FFus = parameters.doubles("FFus");
	double kFus = parameters.doubles("kFus");
	double kDegVen = (1-FFus)/FFus*kFus;
	
	// parameter for time stepping
	double delta_t = 6;// check?
	
	// evaluating Heldt and Laske's ODE models
	pCell->custom_data["VEx"] += (kDislo*VAttlo+kDishi*VAtthi-(kAtthi*Bhi+kAttlo*Blo)*VEx)*delta_t;
	pCell->custom_data["VAtthi"] += (kAtthi*Bhi*VEx-kDishi*VAtthi-kEn*VAtthi)*delta_t;
	pCell->custom_data["VAttlo"] += (kAttlo*Blo*VEx-kDislo*VAttlo-kEn*VAttlo)*delta_t;
	pCell->custom_data["Bhi"] += Btothi-VAtthi;
	pCell->custom_data["Blo"] += Btotlo-VAttlo;
	pCell->custom_data["VEn"] += (kEn*(VAtthi+VAttlo)-(kFus-kDegVen)*VEn)*delta_t;
	
	pCell->phenotype.secretion.uptake_rates[nV_external] = (kAtthi*Bhi+kAttlo*Blo)*Vvoxel;
	pCell->phenotype.secretion.secretion_rates[nV_external] = kDislo*VAttlo+kDishi*VAtthi;
	
	//**************************************************************************
	

	/*
	OLD VIRUS UPTAKE MODEL 
	
	if( rho_virus<1e-6)
	{
		pCell->phenotype.secretion.uptake_rates[nV_external] = 0;
	}
	else if(pCell->custom_data["antiviral_state"]>0.5)
	{
		pCell->phenotype.secretion.uptake_rates[nV_external] = 0;	
	}
	else if(rho_virus>rho_max)
	{
		pCell->phenotype.secretion.uptake_rates[nV_external] = u_Evirus_IFN*(rho_max*rho_max/(m_i/Vvoxel+m_half/Vvoxel))*1/rho_virus;
	}
	else
	{ 
		//pCell->phenotype.secretion.uptake_rates[nV_external] = u_Evirus_IFN*(rho_virus/(rho_half+rho_virus))*(m_half/(m_i+m_half));
		pCell->phenotype.secretion.uptake_rates[nV_external] = u_Evirus_IFN*(rho_virus/(m_i/Vvoxel+m_half/Vvoxel));

	}
	*/
	//****************************************************************
		
		
		
	//if( rho_virus >100)
	//{std::cout<<"Extracellular virus: "<<rho_virus<<" rho_half "<<rho_half<<" rho_virus "<<rho_virus<<" m_half "<<m_half<<" m_i "<<m_i<<std::endl;}

	//if(m_i>7000)
	//{std::cout<<" "<<(rho_virus/(rho_half+rho_virus))<<" "<<(m_half/(m_i+m_half))<<", Intracellular virus: "<<m_i<<std::endl;}

	//if( (rho_virus/(rho_half+rho_virus))*(m_half/(m_i+m_half)) >1 )
	//{std::cout<<"Uptake rate: "<<(rho_virus/(rho_half+rho_virus))*(m_half/(m_i+m_half))<<std::endl;}

	return;
}

void simple_receptor_dynamics_main_model( double dt )
{
	//std::cout<<"simple receptor dynamics main model"<<std::endl;
	#pragma omp parallel for 
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( pC->phenotype.death.dead == false )
		{ 
			simple_receptor_dynamics_model( pC, pC->phenotype , dt ); 
		}
		//std::cout<<"here 3"<<std::endl;
	}
	//std::cout<<"here 2"<<std::endl;
	
	return; 
}
