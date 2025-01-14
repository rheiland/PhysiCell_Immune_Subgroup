#include "./internal_viral_response.h" 

using namespace PhysiCell; 

std::string internal_virus_response_version = "0.4.0"; 

Submodel_Information internal_virus_response_model_info; 

void simple_internal_virus_response_model_setup( void )
{
	// set up the model 
		// set version info 
	internal_virus_response_model_info.name = "internal viral response"; 
	internal_virus_response_model_info.version = internal_virus_response_version; 
		// set functions 
	internal_virus_response_model_info.main_function = NULL; 
	internal_virus_response_model_info.phenotype_function = simple_internal_virus_response_model; 
	internal_virus_response_model_info.mechanics_function = NULL; 
	
	// what microenvironment variables do you need 
	internal_virus_response_model_info.microenvironment_variables.push_back( "virion" ); 
	internal_virus_response_model_info.microenvironment_variables.push_back( "VTEST" ); 		
	internal_virus_response_model_info.microenvironment_variables.push_back( "interferon 1" ); 	
	internal_virus_response_model_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" ); 	
	internal_virus_response_model_info.microenvironment_variables.push_back( "debris" ); 	
	internal_virus_response_model_info.microenvironment_variables.push_back( "chemokine" ); 			
		// what microenvironment variables do you expect? 
		// what custom data do I need? 
	internal_virus_response_model_info.cell_variables.push_back( "r_max" ); 
	internal_virus_response_model_info.cell_variables.push_back( "max_apoptosis_half_max" ); 
	internal_virus_response_model_info.cell_variables.push_back( "apoptosis_hill_power" ); 	
	internal_virus_response_model_info.cell_variables.push_back( "infected_cell_chemokine_secretion_activated" );
	internal_virus_response_model_info.cell_variables.push_back( "infected_cell_chemokine_secretion_rate" );
	internal_virus_response_model_info.cell_variables.push_back( "activated_cytokine_secretion_rate" );
	
	
	internal_virus_response_model_info.cell_variables.push_back( "VEn" ); 	
	internal_virus_response_model_info.cell_variables.push_back( "Vnuc" ); 
	internal_virus_response_model_info.cell_variables.push_back( "VRel" ); 
	
		// register the submodel  
	internal_virus_response_model_info.register_model();	
		// set functions for the corresponding cell definition 
		
//	pCD = find_cell_definition( "lung epithelium" ); 
//	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
//	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	
	return; 
}

void simple_internal_virus_response_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	
		double Vnuc = pCell->custom_data["Vnuc"];
	if(pCell->custom_data["Vnuc"]<1e-6)
	{return;}
	
	static Cell_Definition* pCD = find_cell_definition( "lung epithelium" ); 
	
	// bookkeeping -- find microenvironment variables we need
	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int vtest_external = microenvironment.find_density_index( "VTEST" ); 
	
	static int chemokine_index = microenvironment.find_density_index( "chemokine" );
	static int IFN_index = microenvironment.find_density_index( "interferon 1" );
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
		
	simple_viral_secretion_model( pCell, phenotype, dt );
		
	// over the life time of infection there is a small probability of death
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "apoptosis" );
	phenotype.death.rates[apoptosis_model_index] = parameters.doubles("r_max");//base_death_rate + additional_death_rate; 
		
				
	if( Vnuc > parameters.doubles("Infection_detection_threshold")/Vvoxel - 1e-16 ) 
	{
		pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 1.0; 
	}

	if( pCell->custom_data["infected_cell_chemokine_secretion_activated"] > 0.5 && phenotype.death.dead == false )
	{
		//phenotype.secretion.secretion_rates[chemokine_index] = pCell->custom_data[ "infected_cell_chemokine_secretion_rate" ];//rate; 
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0;
	
		// if cell is infected, start secretion
		phenotype.secretion.secretion_rates[IFN_index]=parameters.doubles("IFN_secretion_rate");

		// (Adrianne) adding pro-inflammatory cytokine secretion by infected cells
		if(pCell->custom_data["antiviral_state"]<0.5)
		{
			double rate = Vnuc*Vvoxel;
			rate /= pCell->custom_data["max_apoptosis_half_max"];
			if(rate>1.0)
			{rate = 1;}
			rate *= pCell->custom_data["infected_cell_chemokine_secretion_rate"];
			
			pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = pCell->custom_data["activated_cytokine_secretion_rate"];
			pCell->phenotype.secretion.secretion_rates[chemokine_index] = rate;
		}
		
		//phenotype.secretion.secretion_rates[vtest_external]  = 1;
	}
	
	if(pCell->custom_data["antiviral_state"]>0.5)
	{
			pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
			pCell->phenotype.secretion.secretion_rates[chemokine_index] = 0;
			pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 0.5; 
	}
		
	return; 	
}

void simple_viral_secretion_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	//std::cout<<pCell->phenotype.secretion.secretion_rates[nV_external]<<std::endl;
	
		static int nV_external = microenvironment.find_density_index( "virion" ); 
		static int vtest_external = microenvironment.find_density_index( "VTEST" ); 
		static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
			
		double Vvoxel = microenvironment.mesh.voxels[1].volume;
		
		if(pCell->phenotype.molecular.internalized_total_substrates[vtest_external]*Vvoxel>8e3 && PhysiCell_globals.current_time>pCell->custom_data["eclipse_time"])
		{
			pCell->phenotype.secretion.secretion_rates[vtest_external] = parameters.doubles("kRel");
				//std::cout<<pCell->phenotype.molecular.internalized_total_substrates[vtest_external]<<std::endl;
		}
		else
		{pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;}	
	
		if(pCell->custom_data["antiviral_state"]>0.5) // cell is in an antiviral state
		{			
			pCell->phenotype.secretion.secretion_rates[vtest_external] = 0;
			pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
			pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
			pCell->custom_data["Vnuc"] = 0;
		}
		
		/*
		double VEx = pCell->nearest_density_vector()[vtest_external]*Vvoxel;//pCell->custom_data["VEx"];
		if(VEx<0.9)
		{VEx = 0;}
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
		
		pCell->custom_data["VRel"] += (kRel*Vnuc)*dt;
		//pCell->phenotype.secretion.secretion_rates[nV_external] = kDislo*VAttlo+kDishi*VAtthi+kRel*Vnuc/Vvoxel; // does this need a dt or virions in the first two terms?
		
		pCell->phenotype.secretion.secretion_rates[vtest_external] = (kDislo*VAttlo+kDishi*VAtthi+kRel*Vnuc)/Vvoxel; // does this need a dt or virions in the first two terms?
		//pCell->phenotype.secretion.secretion_rates[vtest_external] = kRel*Vnuc/Vvoxel; // does this need a dt or virions in the first two terms?
		
		*/
	return;
}




