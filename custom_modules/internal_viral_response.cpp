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
	
		// what microenvironment variables do you expect? 
		// what custom data do I need? 
	internal_virus_response_model_info.cell_variables.push_back( "max_infected_apoptosis_rate" ); 
	internal_virus_response_model_info.cell_variables.push_back( "max_apoptosis_half_max" ); 
	internal_virus_response_model_info.cell_variables.push_back( "apoptosis_hill_power" ); 	
	
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
	static Cell_Definition* pCD = find_cell_definition( "lung epithelium" ); 
	
	// bookkeeping -- find microenvironment variables we need

	static int nV_external = microenvironment.find_density_index( "virion" ); 
	static int nA_external = microenvironment.find_density_index( "assembled virion" ); 
	static int chemokine_index = microenvironment.find_density_index( "chemokine" );
	static int IFN_index = microenvironment.find_density_index( "interferon 1" );
	
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
	
	static int nV_internal = pCell->custom_data.find_variable_index( "virion" ); 
	 
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "apoptosis" );
	
	double R = pCell->phenotype.molecular.internalized_total_substrates[nV_external];
	
	// base death rate (from cell line)
	double base_death_rate = 
		pCD->phenotype.death.rates[apoptosis_model_index]; 
	
	// additional death rate from infectoin  
	double additional_death_rate = pCell->custom_data["max_infected_apoptosis_rate"] ; 
	
	double v = pCell->phenotype.molecular.internalized_total_substrates[nV_external] / 
		pCell->custom_data["max_apoptosis_half_max"] ; 
	v = pow( v, pCell->custom_data["apoptosis_hill_power"] ); 
	
	double effect = v / (1.0+v); 
	additional_death_rate *= effect; 
	phenotype.death.rates[apoptosis_model_index] = base_death_rate + additional_death_rate; 
	
	// if we're infected, secrete a chemokine for the immune model
	double AV = pCell->phenotype.molecular.internalized_total_substrates[nV_external];  
	
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
		
	
	if( R >= parameters.doubles("Infection_detection_threshold")/Vvoxel - 1e-16 ) 
	{
		pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 1.0; 
	}

	if( pCell->custom_data["infected_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
	{
		double rate = AV; 
		rate /= pCell->custom_data["max_apoptosis_half_max"];
		if( rate > 1.0 )
		{ rate = 1.0; }
		rate *= pCell->custom_data[ "infected_cell_chemokine_secretion_rate" ];

		phenotype.secretion.secretion_rates[chemokine_index] = pCell->custom_data[ "infected_cell_chemokine_secretion_rate" ];//rate; 
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0;
	
		// if cell is infected, start secretion
		phenotype.secretion.secretion_rates[IFN_index]=parameters.doubles("IFN_secretion_rate");

		// (Adrianne) adding pro-inflammatory cytokine secretion by infected cells
		phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = pCell->custom_data["activated_cytokine_secretion_rate"];
	}
		
	return; 	
}