#include "./epithelium_submodel.h" 

using namespace PhysiCell; 

std::string epithelium_submodel_version = "0.4.0"; 

Submodel_Information epithelium_submodel_info; 

void epithelium_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
	// elastic adhesions 
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );
	
	return; 
}

void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	
	// receptor dynamics 
	// requires faster time scale - done in main function 
	
	// viral dynamics model 
	internal_viral_dynamics_info.phenotype_function(pCell,phenotype,dt); 
	// internal_virus_model(pCell,phenotype,dt);
	
	// viral response model 
	internal_virus_response_model_info.phenotype_function(pCell,phenotype,dt); 
	// internal_virus_response_model(pCell,phenotype,dt);	
	
	// T-cell based death
	TCell_induced_apoptosis(pCell, phenotype, dt ); 
	
	// (Adrianne) ROS induced cell death model
	ROS_induced_apoptosis(pCell, phenotype, dt);
	
	// (Adrianne) cell proliferation model
	Cell_proliferation(pCell, phenotype, dt);
		
	// if I am dead, remove all adhesions 
	static int apoptosis_index = phenotype.death.find_death_model_index( "apoptosis" ); 
	if( phenotype.death.dead == true )
	{
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		
		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
	}

	static int IFN_index = microenvironment.find_density_index( "interferon 1" );
	double IFN_internal = pCell->nearest_density_vector()[IFN_index];
	double IC_50_IFN = parameters.doubles("IC_50_IFN");
	double Vvoxel = microenvironment.mesh.voxels[1].volume;
	
	double IFN_prob = IFN_internal/(IC_50_IFN+IFN_internal);
	
	double prob_prob = UniformRandom();
	if( prob_prob<IFN_prob)
	{
		// cell enters antiviral state
		//std::cout<<"Cell enters antiviral state "<<prob_prob<<" "<<IFN_prob<<std::endl;
		pCell->custom_data["antiviral_state"] = 1;
	}
	
/*
	// cell secretion belongs in viral response 
	
	// if I am dead, make sure to still secrete the chemokine 
	static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
	static int nP = pCell->custom_data.find_variable_index( "viral_protein"); 
	double P = pCell->custom_data[nP];
	
	static int nAV = pCell->custom_data.find_variable_index( "assembled_virion" ); 
	double AV = pCell->custom_data[nAV]; 

	static int nR = pCell->custom_data.find_variable_index( "viral_RNA");
	double R = pCell->custom_data[nR];
	
	if( R >= 1.00 - 1e-16 ) 
	{
		pCell->custom_data["infected_cell_chemokine_secretion_activated"] = 1.0; 
	}

	if( pCell->custom_data["infected_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
	{
		double rate = AV; // P; 
		rate /= pCell->custom_data["max_apoptosis_half_max"];
		if( rate > 1.0 )
		{ rate = 1.0; }
		rate *= pCell->custom_data[ "infected_cell_chemokine_secretion_rate" ];

		phenotype.secretion.secretion_rates[chemokine_index] = rate; 
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0; 
	}
*/	
	
	// if I am dead, don't bother executing this function again 
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}

void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	
	pCell->is_movable = false; 
	
	// if I'm dead, don't bother 
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions, 
		// since those are part of mechanics. 
		// remove_all_adhesions( pCell ); 
		
		// Let's just fully disable now. 
		pCell->functions.custom_cell_rule = NULL; 
		pCell->functions.contact_function = NULL; 

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		return; 
	}	
	
	// this is now part of contact_function 
	/*
	// if I'm adhered to something ... 
	if( pCell->state.neighbors.size() > 0 )
	{
		// add the elastic forces 
		extra_elastic_attachment_mechanics( pCell, phenotype, dt );
	}
	*/
	return; 
}

void epithelium_submodel_setup( void )
{
	Cell_Definition* pCD;
	
	// set up any submodels you need 
	// viral replication 
	
	// receptor trafficking 
	simple_receptor_dynamics_model_setup(); // done 
	// viral replication 
	simple_internal_virus_model_setup();	
	// single-cell response  XXXXXX THIS NEEDS TO BE EDITED STILL XXXXXXXXXXXXXXXXXXXXXXX
	simple_internal_virus_response_model_setup(); 
 	
	// set up epithelial cells
		// set version info 
	epithelium_submodel_info.name = "epithelium model"; 
	epithelium_submodel_info.version = epithelium_submodel_version; 
		// set functions 
	epithelium_submodel_info.main_function = NULL; 
	epithelium_submodel_info.phenotype_function = epithelium_phenotype; 
	epithelium_submodel_info.mechanics_function = epithelium_mechanics; 
	
		// what microenvironment variables do you expect? 
	epithelium_submodel_info.microenvironment_variables.push_back( "virion" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "interferon 1" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "chemokine" ); 
		// what custom data do I need? 
	//epithelium_submodel_info.cell_variables.push_back( "something" ); 
		// register the submodel  
	epithelium_submodel_info.register_model();	
		// set functions for the corresponding cell definition 
	pCD = find_cell_definition( "lung epithelium" ); 
	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	pCD->functions.contact_function = epithelium_contact_function; 
	
	return;
}

void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" ); 
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
	static int virion_index = microenvironment.find_density_index("virion");
	
	
	if( pCell->custom_data["TCell_contact_time"] > pCell->custom_data["TCell_contact_death_threshold"] )
	{
		// make sure to get rid of all adhesions! 
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		
		#pragma omp critical
		{
		std::cout << "\t\t\t\t" << pCell << " (of type " << pCell->type_name <<  ") died from T cell contact" << std::endl; 
		}
		
		// induce death 
		pCell->start_death( apoptosis_index ); 
		
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0; 
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		pCell->phenotype.molecular.fraction_released_at_death[virion_index] = parameters.doubles("virus_fraction_released_after_apoptosis");
		
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}


void ROS_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" ); 
	static int ROS_index = microenvironment.find_density_index( "ROS" ); 
	double ROS_amount = pCell->nearest_density_vector()[ROS_index];
	static int debris_index = microenvironment.find_density_index( "debris" ); 
	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
	
	double epsilon_ROS = parameters.doubles("epsilon_ROS");
	
	double prob_apoptosis = ROS_amount/(ROS_amount+epsilon_ROS);
	
	if( UniformRandom() < prob_apoptosis )
	{
		// make sure to get rid of all adhesions! 
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		
		#pragma omp critical
		{
		std::cout << "\t\t\t\t" << pCell << " (of type " << pCell->type_name <<  ") died from ROS" << std::endl; 
		}
		
		// induce death 
		pCell->start_death( apoptosis_index ); 
		
		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0; 
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
		
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}

void Cell_proliferation( Cell* pCell, Phenotype& phenotype, double dt )
{
		if(pCell->custom_data["cell_attachment_lifetime"]>50)
	{std::cout<<"Cell attachment: "<<pCell->custom_data["cell_attachment_lifetime"]<<std::endl;
		pCell->die();}
			
	/*int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
	
	double c_j_ccr = 10;
	double cell_radius = 8.4;
	double A_cell = 4*3.141*cell_radius*cell_radius;
	double A_frag = 800*800;
	double K = 2793;
	double pressure_threshold = 6.98;//6*c_j_ccr/A_cell*(1-1/(2*cell_radius)*sqrt(2*A_frag/(sqrt(3)*K)))/(0.02729*c_j_ccr/A_cell);
	
	// first - cells on boundary don't proliferate!
	std::vector<double> cells_position(3);
    cells_position = pCell->position;


	if(cells_position[0]>default_microenvironment_options.X_range[1]-pCell->phenotype.geometry.radius*3
	|| cells_position[0]<default_microenvironment_options.X_range[0]+pCell->phenotype.geometry.radius*3
	|| cells_position[1]>default_microenvironment_options.Y_range[1]-pCell->phenotype.geometry.radius*3
	|| cells_position[1]<default_microenvironment_options.Y_range[0]+pCell->phenotype.geometry.radius*3)
	{	
		pCell->phenotype.cycle.data.transition_rate( G0G1_index , S_index ) = 0;
		return;
	}
			
	// if not boundary cell, and the simp press is less than threshold then proliferate
	double simp_press = pCell->state.simple_pressure;//calculate_simple_pressure_again(pCell);
	if(simp_press<0.53 )//parameters.doubles("pressure_threshold"))
	{
		if(simp_press>1e-6)
		{//std::cout<<simp_press<<std::endl;
		pCell->phenotype.cycle.data.transition_rate( G0G1_index , S_index ) = parameters.doubles("epithelial_cell_proliferation_rate");
		pCell->is_movable=true;
		pCell->phenotype.motility.migration_speed = 0.05;
		}
	}
	else
	{
		pCell->phenotype.cycle.data.transition_rate( G0G1_index , S_index ) = 0.0;
		//pCell->is_movable=false;
		//pCell->phenotype.motility.migration_speed = 0.01;
	}
	*/
	return;

}

double calculate_simple_pressure_again(Cell* pCell)
{
	
	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2''
	
	std::vector<Cell*> neighbors = pCell->cells_in_my_container(); 

	// if the only cell is pCell, return simp_press of 0 (as in there is no pressure so cell can proliferate)
	if( neighbors.size() < 2 )
	{return 0.0; } 

	double simp_press = 0; 
	int n = 0; 
	Cell* pContactCell = neighbors[n]; 
	while( n < neighbors.size() )
	{
		pContactCell = neighbors[n]; 
		if(pContactCell->type == pCell->type)
		{	double distance = 1-sqrt((pContactCell->position[0]-pCell->position[0])*(pContactCell->position[0]-pCell->position[0])+(pContactCell->position[1]-pCell->position[1])*(pContactCell->position[1]-pCell->position[1]));
			double distance_norm_to_radius = distance/(2*pCell->phenotype.geometry.radius);
			double displacement = 1-distance_norm_to_radius;
			simp_press += displacement*displacement;
		}
		else		
		{simp_press +=0;}
		n++;
	}
	
	pCell->custom_data["displacement_stor"] = simp_press;
		
	return simp_press;

}