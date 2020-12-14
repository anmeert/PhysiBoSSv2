
#include "./custom.h"

/**
 *	\main drug_AGS custom
 *	\brief Custom module file for prostate example
 * 
 *	\details Modules needed for the prostate example. This custom module can be used to study the inhibition of AGS cell lines with AKT, beta-catenin and TAK inhibitors.
 *
 *
 *	\date 19/10/2020
 *	\author Arnau Montagud, BSC-CNS, with code previously developed by Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

// declare cell definitions here 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.functions.set_orientation = NULL;

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	initialize_cell_definitions_from_pugixml();

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	
	// set molecular properties 
	int erki_substrate_index = microenvironment.find_density_index( "ERKi" );
	cell_defaults.phenotype.molecular.fraction_released_at_death[erki_substrate_index] = 0.0;

	int myc_maxi_substrate_index = microenvironment.find_density_index( "MYC_MAXi" );
	cell_defaults.phenotype.molecular.fraction_released_at_death[myc_maxi_substrate_index] = 0.0;


	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 

	return; 
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	

	// initialize BioFVM 
	initialize_microenvironment(); 	
	
	return; 
}
void update_custom_variables( Cell* pCell )
{
	// static int tnf_index = microenvironment.find_density_index( "tnf" ); 
	// static int index_tnf_concentration = pCell->custom_data.find_variable_index("tnf_concentration");
	// static int index_tnf_node = pCell->custom_data.find_variable_index("tnf_node");
	// static int index_fadd_node = pCell->custom_data.find_variable_index("fadd_node");
	// pCell->custom_data.variables.at(index_tnf_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[tnf_index];
	// pCell->custom_data.variables.at(index_tnf_node).value = pCell->boolean_network.get_node_value("TNF");
	// pCell->custom_data.variables.at(index_fadd_node).value = pCell->boolean_network.get_node_value("FADD");
	
	static int erki_index = microenvironment.find_density_index( "ERKi" ); 
	static int index_erki_concentration = pCell->custom_data.find_variable_index("erki_concentration");
	static int index_erki_node = pCell->custom_data.find_variable_index("erki_node");
	pCell->custom_data.variables.at(index_erki_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[erki_index];
	//pCell->custom_data.variables.at(index_erki_node).value = pCell->boolean_network.get_node_value("anti_ERK");

	static int myc_maxi_index = microenvironment.find_density_index( "MYC_MAXi" ); 
	static int index_myc_maxi_concentration = pCell->custom_data.find_variable_index("myc_maxi_concentration");
	static int index_myc_maxi_node = pCell->custom_data.find_variable_index("myc_maxi_node");
	pCell->custom_data.variables.at(index_myc_maxi_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[myc_maxi_index];
	//pCell->custom_data.variables.at(index_myc_maxi_node).value = pCell->boolean_network.get_node_value("anti_MYC_MAX");

	// static int akti_index = microenvironment.find_density_index( "AKTi" ); 
	// static int index_akti_concentration = pCell->custom_data.find_variable_index("akti_concentration");
	// static int index_akti_node = pCell->custom_data.find_variable_index("akti_node");
	// pCell->custom_data.variables.at(index_akti_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[akti_index];
	// pCell->custom_data.variables.at(index_akti_node).value = pCell->boolean_network.get_node_value("anti_AKT");

	// static int taki_index = microenvironment.find_density_index( "TAKi" ); 
	// static int index_taki_concentration = pCell->custom_data.find_variable_index("taki_concentration");
	// static int index_taki_node = pCell->custom_data.find_variable_index("taki_node");
	// pCell->custom_data.variables.at(index_taki_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[taki_index];
	// pCell->custom_data.variables.at(index_taki_node).value = pCell->boolean_network.get_node_value("anti_TAK1");

	// static int bcati_index = microenvironment.find_density_index( "BCATi" ); 
	// static int index_bcati_concentration = pCell->custom_data.find_variable_index("bcati_concentration");
	// static int index_bcati_node = pCell->custom_data.find_variable_index("bcati_node");
	// pCell->custom_data.variables.at(index_bcati_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[bcati_index];
	// pCell->custom_data.variables.at(index_bcati_node).value = pCell->boolean_network.get_node_value("anti_betacatenin");


}

void setup_tissue( void )
{
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	std::string bnd_file = parameters.strings("bnd_file");
	std::string cfg_file = parameters.strings("cfg_file");
	BooleanNetwork prostate_network;
	double maboss_time_step = parameters.doubles("maboss_time_step");
	prostate_network.initialize_boolean_network(bnd_file, cfg_file, maboss_time_step);

	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		pC = create_cell(get_cell_definition("default")); 
		pC->assign_position( x, y, z );
		// pC->set_total_volume(sphere_volume_from_radius(radius));
		
		// pC->phenotype.cycle.data.current_phase_index = phase;
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		
		pC->boolean_network = prostate_network;
		pC->boolean_network.restart_nodes();
		static int index_next_physiboss_run = pC->custom_data.find_variable_index("next_physiboss_run");
		pC->custom_data.variables.at(index_next_physiboss_run).value = pC->boolean_network.get_time_to_update();
		update_custom_variables(pC);
	}

	return; 
}

// custom cell phenotype function to run PhysiBoSS when is needed
void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	if ( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model) 
	{
		std::cout << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) << "\n" << std::endl;
	}
	// if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
	// {
	// 	std::cout << pCell->phenotype.cycle.current_phase().name << " 0,1: " << pCell->phenotype.cycle.data.transition_rate(0, 1) << " / " << pCell->phenotype.cycle.current_phase().name << " 1,0: " << pCell->phenotype.cycle.data.transition_rate(1, 0)  << "\n" << std::endl;
	// }

	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
	static int index_next_physiboss_run = pCell->custom_data.find_variable_index("next_physiboss_run");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	if (PhysiCell_globals.current_time >= pCell->custom_data.variables.at(index_next_physiboss_run).value)
	{
		set_input_nodes(pCell);

		pCell->boolean_network.run_maboss();
		// Get noisy step size
		double next_run_in = pCell->boolean_network.get_time_to_update();
		pCell->custom_data.variables.at(index_next_physiboss_run).value = PhysiCell_globals.current_time + next_run_in;
		
		update_custom_variables(pCell);

		from_nodes_to_cell(pCell, phenotype, dt);
	}
}

std::vector<std::string> prolif_apoptosis_coloring( Cell* pCell )
{
	std::vector<std::string> output;
	// color apoptotic cells red and proliferating cells green
	if (pCell->boolean_network.get_node_value("Apoptosis"))
	{
		//color red
		output = {"crimson", "black","darkred", "darkred"};
	}
	else 
	{
		//color green
		output = {"limegreen", "black", "darkgreen", "darkgreen"};
	} 
	return output; 
}

void set_boolean_node (Cell* pCell, std::string node_to_be_set, int index, double threshold) {
	if (index != -1)
		{
			double cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[index];
			if (cell_concentration >= threshold)
				pCell->boolean_network.set_node_value(node_to_be_set, 1);
			else 
			{
				pCell->boolean_network.set_node_value(node_to_be_set, 0);
			}
		}
}


void set_input_nodes(Cell* pCell) {
	
	static int erki_index = microenvironment.find_density_index("ERKi");
	static double erki_threshold = parameters.doubles("erki_threshold");
	static int myc_maxi_index = microenvironment.find_density_index("MYC_MAXi");
	static double myc_maxi_threshold = parameters.doubles("myc_maxi_threshold");

	// In order to do that i first need to create .cfg file that contain the anti-nodes 
	//set_boolean_node("anti_ERK", erki_index, erki_threshold);
	//set_boolean_node("anti_MYC_MAX", myc_maxi_index, myc_maxi_threshold);
}


void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	std::vector<bool>* nodes = pCell->boolean_network.get_nodes();
	int bn_index;

	// For prostate model

	// live model 
	// translate apoptosis, proliferation and invasion values into agent-based model

	if( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
	{
		int start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		int end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		int apoptosis_index = phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
		
		if(pCell->boolean_network.get_node_value("Apoptosis"))
		{
			pCell->start_death(apoptosis_index);
			return;
		}

		if( pCell->boolean_network.get_node_value("EMT"))
		{
			// pCell->evolve_coef_sigmoid( 
			// 	pCell->boolean_network.get_node_value("EMT"), phenotype.mechanics.cell_cell_adhesion_strength, dt 
			// );

			//phenotype.mechanics.cell_cell_adhesion_strength = PhysiCell::parameters.doubles("homotypic_adhesion_max") * theta 
		}

		if(pCell->boolean_network.get_node_value("Proliferation"))
		{
			double multiplier = 2.0;
			phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) *= multiplier;
		}

		std::cout << "Migration: " << pCell->boolean_network.get_node_value("Migration") << std::endl;
		
		if(pCell->boolean_network.get_node_value("Migration"))
		{ 
			pCell->phenotype.motility.is_motile = true;
		 	pCell->phenotype.motility.migration_speed = PhysiCell::parameters.doubles("migration_speed");
			pCell->phenotype.motility.migration_bias = PhysiCell::parameters.doubles("migration_bias");
			pCell->phenotype.motility.persistence_time = PhysiCell::parameters.doubles("persistence");

			// pCell->evolve_coef(pCell->boolean_network.get_node_value("Migration"),	phenotype.motility.migration_speed, dt 
			// );
			// pCell->phenotype.motility.migration_speed = PhysiCell::parameters.doubles("max_motility_speed") * migration_coeff
		}

		if(pCell->boolean_network.get_node_value("Invasion"))
		{
			// nothing happens for now 
		}

	}


	// if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
	// {
	// 	std::cout << pCell->phenotype.cycle.current_phase().name << " 0,1: " << pCell->phenotype.cycle.data.transition_rate(0, 1) <<  ". " << prosurvival_value << ", " << antisurvival_value << " / " << pCell->phenotype.cycle.current_phase().name << " 1,0: " << pCell->phenotype.cycle.data.transition_rate(1, 0) <<  ". " << prosurvival_value << ", " << antisurvival_value << "\n" << std::endl;

	// }

	// int tnf_substrate_index = microenvironment.find_density_index( "tnf" );
	// static double tnf_secretion = parameters.doubles("tnf_secretion_rate");
	// double tnf_secretion_rate = 0;
	// // produce some TNF
	// if ( pCell->boolean_network.get_node_value( "NFkB" ) )
	// {
	// 	tnf_secretion_rate = (tnf_secretion / microenvironment.voxels(pCell->get_current_voxel_index()).volume);

	// }
	// pCell->phenotype.secretion.secretion_rates[tnf_substrate_index] = tnf_secretion_rate;
	pCell->set_internal_uptake_constants(dt);
}

// ***********************************************************
// * NOTE: Funtion replicated from PhysiBoSS, but not used   *
// *       as we use a live cycle model instead a Ki67 model *
// ***********************************************************
void do_proliferation( Cell* pCell, Phenotype& phenotype, double dt )
{
		//if cells in basic_Ki67_cycle_model (code 1)
	// TODO adaptar rate conforme la resta Prosurv - Antisurv
		// If cells is in G0 (quiescent) switch to pre-mitotic phase
		// if ( pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::Ki67_negative ) {
		// if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
		

	// if cells in live_cells_cycle_model (code 5)
	// TODO adaptar rate conforme la resta Prosurv - Antisurv
		// if( phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
		// if ( pCell->phenotype.cycle.current_phase_index() == PhysiCell_constants::live) {
			// pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt);
		// }
}

// ***********************************************************
// * NOTE: Funtion to read init files created with PhysiBoSS *
// ***********************************************************
std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header) 
{ 
	// File pointer 
	std::fstream fin; 
	std::vector<init_record> result;

	// Open an existing file 
	fin.open(filename, std::ios::in); 

	// Read the Data from the file 
	// as String Vector 
	std::vector<std::string> row; 
	std::string line, word;

	if(header)
		getline(fin, line);

	do 
	{
		row.clear(); 

		// read an entire row and 
		// store it in a string variable 'line' 
		getline(fin, line);

		// used for breaking words 
		std::stringstream s(line); 

		// read every column data of a row and 
		// store it in a string variable, 'word' 
		while (getline(s, word, delimiter)) { 

			// add all the column data 
			// of a row to a vector 
			row.push_back(word); 
		}

		init_record record;
		record.x = std::stof(row[2]);
		record.y = std::stof(row[3]);
		record.z = std::stof(row[4]);
		record.radius = std::stof(row[5]);
		record.phase = std::stoi(row[13]);
		record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());
	
	return result;
}