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
	
	
	// what microenvironment variables do you need 
	internal_viral_dynamics_info.microenvironment_variables.push_back( "virion" ); 	
	internal_viral_dynamics_info.microenvironment_variables.push_back( "interferon 1" ); 		
	// what microenvironment variables do I need? 

	// what custom data do I need? 
	internal_viral_dynamics_info.cell_variables.push_back( "VEn" ); 		
	internal_viral_dynamics_info.cell_variables.push_back( "Vnuc" ); 
		
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
					
	simple_intracellular_replication_model(pCell, phenotype, dt);
			
	return;
}
void simple_intracellular_replication_model(  Cell* pCell, Phenotype& phenotype, double dt )
{	

	static int vtest_external = microenvironment.find_density_index( "VTEST" ); 
	
	static int proinflammatory_cytokine_index = microenvironment.find_density_index( "pro-inflammatory cytokine");
	
	double VEn = pCell->custom_data["VEn"];
	
	double Vconc = pCell->phenotype.molecular.internalized_total_substrates[vtest_external];
		
		double Vvoxel = microenvironment.mesh.voxels[1].volume;
	//static double kFus = parameters.doubles("kFus");
	static double gamnuc = parameters.doubles("gamnuc");
	static double KVnuc = parameters.doubles("KVnuc")/Vvoxel;
	//static double kRel = parameters.doubles("kRel");
	
	//std::cout<<VEn<<" "<<Vnuc<<std::endl;
	if(pCell->custom_data["antiviral_state"]<0.5&&pCell->phenotype.molecular.internalized_total_substrates[vtest_external]*Vvoxel>10) // cell isn't in an antiviral state
	{	
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] += gamnuc*Vconc*(1-Vconc/KVnuc)*dt;
		if(pCell->custom_data["eclipse_time"]<1)
		{pCell->custom_data["eclipse_time"] = PhysiCell_globals.current_time+12*60;}
	}
	else if(pCell->custom_data["antiviral_state"]>0.5)
	{
		pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;
		phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
	}
	else if(pCell->phenotype.molecular.internalized_total_substrates[vtest_external]*Vvoxel<11)
	{
			double prob_infection_recognition = 0.3; // probability at low MOI a cell realises it's infected
			if(UniformRandom()<prob_infection_recognition)
			{pCell->phenotype.molecular.internalized_total_substrates[vtest_external] = 0;}
			
	}
	pCell->custom_data["Vnuc"] = pCell->phenotype.molecular.internalized_total_substrates[vtest_external];	
	
	return;	
}
