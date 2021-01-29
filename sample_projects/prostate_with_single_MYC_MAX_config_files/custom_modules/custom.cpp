
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
	
	// int erki_substrate_index = microenvironment.find_density_index( "ERKi" );
	// cell_defaults.phenotype.molecular.fraction_released_at_death[erki_substrate_index] = 0.0;

	// int myc_maxi_substrate_index = microenvironment.find_density_index( "MYC_MAXi" );
	// cell_defaults.phenotype.molecular.fraction_released_at_death[myc_maxi_substrate_index] = 0.0;

	// set molecular properties for each drug defined in PhysiCell_settings.xml
	// starts with the second density because the first density is oxygen
	for (int i = 1; i < microenvironment.number_of_densities(); i++) 
	{
		cell_defaults.phenotype.molecular.fraction_released_at_death[i] = 0.0;
	}

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
	// first density is oxygen - shouldn't be changed: index from 1
	for (int i = 1; i < microenvironment.number_of_densities(); i++) 
	{
		std::string drug_name = microenvironment.density_names[i];
		int drug_index = microenvironment.find_density_index(drug_name);
		int index_drug_conc = pCell->custom_data.find_variable_index(drug_name + "_concentration");
		int index_drug_node = pCell->custom_data.find_variable_index(drug_name + "_node");
		pCell->custom_data.variables.at(index_drug_conc).value = pCell->phenotype.molecular.internalized_total_substrates[drug_index];
		pCell->custom_data.variables.at(index_drug_node).value = pCell->boolean_network.get_node_value("anti_" + drug_name);

		
	}
	// static int erki_index = microenvironment.find_density_index( "ERKi" ); 
	// static int index_erki_concentration = pCell->custom_data.find_variable_index("erki_concentration");
	// static int index_erki_node = pCell->custom_data.find_variable_index("erki_node");
	// pCell->custom_data.variables.at(index_erki_concentration).value = pCell->phenotype.molecular.internalized_total_substrates[erki_index];
	//pCell->custom_data.variables.at(index_erki_node).value = pCell->boolean_network.get_node_value("anti_ERK");
}

void setup_tissue( void )
{
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	std::string bnd_file = PhysiCell::parameters.strings("bnd_file");
	std::string cfg_file = PhysiCell::parameters.strings("cfg_file");
	BooleanNetwork prostate_network;
	double maboss_time_step = PhysiCell::parameters.doubles("maboss_time_step");
	prostate_network.initialize_boolean_network(bnd_file, cfg_file, maboss_time_step);

	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		double random_num_1 = (double) rand()/RAND_MAX;
		double random_num_2 = (double) rand()/RAND_MAX;

		if (PhysiCell::parameters.ints("simulation_mode") == 0)
		{
			// single inhibition - just one drug is present 
			if (random_num_1 < PhysiCell::parameters.doubles("prop_drug_sensitive_" + microenvironment.density_names[1]))
			{
				// cell is sensitive to the drug
				pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_sensitive"));
			}
			else 
			{
				// cell is not sensitive to the drug
				pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_insensitive"));
			}
		}
		else
		{
			// double inhibition - two drugs are present - we have 4 cell strains 
			if (random_num_1 < PhysiCell::parameters.doubles("prop_drug_sensitive_" + microenvironment.density_names[1]))
			{
				if (random_num_2 < PhysiCell::parameters.doubles("prop_drug_sensitive_" + microenvironment.density_names[2]))
				{
					// cell is sensitive to both drugs
					pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_sensitive"));
				}
				else 
				{
					// cell is only sensitive to the first drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[2] + "_insensitive"));
				}
			}
			else
			{
				if (random_num_2 < PhysiCell::parameters.doubles("prop_drug_sensitive_" + microenvironment.density_names[2]))
				{
					// cell is only sensitive to the second drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[2] + "_sensitive"));
				}
				else
				{
					// cell is sensitive to no drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_insensitive"));
				}
				
			}
			
		}
 
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
		//std::cout << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) << "\n" << std::endl;
	}
	// if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
	// {
	// 	std::cout << pCell->phenotype.cycle.current_phase().name << " 0,1: " << pCell->phenotype.cycle.data.transition_rate(0, 1) << " / " << pCell->phenotype.cycle.current_phase().name << " 1,0: " << pCell->phenotype.cycle.data.transition_rate(1, 0)  << "\n" << std::endl;
	// }

	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);

	// update motility state variable
	static int index_motility_state = pCell->custom_data.find_variable_index("motility_state");
	pCell->custom_data.variables.at(index_motility_state).value = int(pCell->phenotype.motility.is_motile);

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
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptosis_death_model)
	{
		//apoptotic cells are colored red
		output = {"crimson", "black","darkred", "darkred"};
	}

	else if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrosis_death_model)
	{
		//necrotic cells are colored brown
		output = {"peru", "black","saddlebrown", "saddlebrown"};
	}

	else if (PhysiCell::parameters.ints("simulation_mode") == 0) 
	{
		std::string drug_name = microenvironment.density_names[1];
		if (pCell->type_name == drug_name + "_sensitive")
		{
			//drug sensitive living cells are colored blue
			output = {"deepskyblue", "black", "darkblue", "darkblue"};
		} 
		else 
		{
			//drug insensitive living cells are colored green
			output = {"limegreen", "black", "darkgreen", "darkgreen"};
		}
	}
	else 
	{
		// // color living cells just in one color 
		// output = {"limegreen", "black", "darkgreen", "darkgreen"};

		// In case we want to color all 4 strains differently:
		// double inhibitions --> 4 cell strains 
		std::string drug1_name = microenvironment.density_names[1];
		std::string drug2_name = microenvironment.density_names[2]; 
		if (pCell->type_name == drug1_name + "_sensitive")
		{
			//cells that are sensitive to both drugs are colored blue
			output = {"deepskyblue", "black", "darkblue", "darkblue"};
		}
		else if (pCell->type_name == drug2_name + "_insensitive")
		{
			// cells that are just sensitive to the first drug 
			output = {"limegreen", "black", "darkgreen", "darkgreen"};

		}
		else if (pCell->type_name == drug2_name + "_sensitive")
		{
			// cells are just sensitive to the second drug
			output = {"gold", "black", "orange", "orange"};
		}
		else 
		{
			// cells aren't sensitive to any drug
			output = {"mediumorchid", "black", "mediumpurple", "mediumpurple"};
		}
	}
	return output;

}

void set_boolean_node (Cell* pCell, std::string node, int index, double threshold) {
	if (index != -1)
		{
			double cell_concentration = pCell->phenotype.molecular.internalized_total_substrates[index];
			if (cell_concentration >= threshold)
			{
				pCell->boolean_network.set_node_value(node, 1);
			}
			else 
			{
				pCell->boolean_network.set_node_value(node, 0);
			}
			
		}
}


void set_input_nodes(Cell* pCell) {

	if (PhysiCell::parameters.ints("simulation_mode") == 0)
	{	
		// single inhibition - just one drug is present 

		std::string drug_name = microenvironment.density_names[1];
		int drug_index = microenvironment.find_density_index(drug_name);
		int drug_threshold = PhysiCell::parameters.doubles( "threshold_" + drug_name);

		if (pCell->type_name == drug_name + "_sensitive")
		{
			// cell is sensitive to the drug -> set boolean node
			set_boolean_node(pCell, "anti_" + drug_name, drug_index, drug_threshold);
		}
	}
	else
	{	
		// double inhibition - two drugs are present - we have 4 cell strains 

		std::string drug1_name = microenvironment.density_names[1];
		int drug1_index = microenvironment.find_density_index(drug1_name);
		int drug1_threshold = PhysiCell::parameters.doubles( "threshold_" + drug1_name);

		std::string drug2_name = microenvironment.density_names[2];
		int drug2_index = microenvironment.find_density_index(drug2_name);
		int drug2_threshold = PhysiCell::parameters.doubles( "threshold_" + drug2_name);
	
		if (pCell->type_name == drug1_name + "_sensitive")
		{	
			// cell is sensitive to both drugs
			set_boolean_node(pCell, "anti_" + drug1_name, drug1_index, drug1_threshold);
			set_boolean_node(pCell, "anti_" + drug2_name, drug2_index, drug2_threshold);
		}
		else if (pCell->type_name == drug2_name + "_insensitive")
		{
			// cell is only sensitive to the first drug
			set_boolean_node(pCell, "anti_" + drug1_name, drug1_index, drug1_threshold);
		}
		else if (pCell->type_name == drug2_name +  "_sensitive")
		{
			// cell is only sensitive to the second drug
			set_boolean_node(pCell, "anti_" + drug2_name, drug2_index, drug2_threshold);
		}

		// else: cell is sensitive to no drug --> no boolean node is set
	}
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
		
		// Update Apoptosis 
		if(pCell->boolean_network.get_node_value("Apoptosis"))
		{
			pCell->start_death(apoptosis_index);
			return;
		}


		// Update Adhesion
		if( pCell->boolean_network.get_node_value("EMT"))
		{
			// reduce cell-cell adhesion 
			// pCell->evolve_coef_sigmoid( 
			// 	pCell->boolean_network.get_node_value("EMT"), phenotype.mechanics.cell_cell_adhesion_strength, dt 
			// );

			//phenotype.mechanics.cell_cell_adhesion_strength = PhysiCell::parameters.doubles("homotypic_adhesion_max") * theta 
		}


		// Update pReference_live_phenotype for proliferation node 
		double max_trans_rate = PhysiCell::parameters.doubles("max_transition_rate");
		double min_trans_rate = PhysiCell::parameters.doubles("min_transition_rate");

		if (pCell->boolean_network.get_node_value("Proliferation")) 
		{
			// multiplier implementation
			//pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) *= 2.5;
			//std::cout << "Rate up! " << std::endl;


			//switch implementation
			pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = PhysiCell::parameters.doubles("high_transition_rate");
		}
		else 
		{
			//multiplier implementation 
			//pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) *= 0.4;
			//std::cout << "Rate down! " << std::endl;


			//switch implementation 
			pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = PhysiCell::parameters.doubles("low_transition_rate");
		}

		 
	
		// if (pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) > max_trans_rate) 
		// {
		// 	 pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = max_trans_rate;
		// }
		// else if (pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) < min_trans_rate)
		// {
		// 	 pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = min_trans_rate;
		// }


		// Update Migration
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
		else 
		{
			pCell->phenotype.motility.is_motile = false;
		}


		// Update Invasion
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
