/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based; 
	
	// only needed for a 2-D simulation: 
	
	/*
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	*/
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 

	// initially no necrosis 
	int live_index = cell_defaults.phenotype.cycle.model().find_phase_index(PhysiCell_constants::live);
	cell_defaults.phenotype.cycle.data.transition_rate(live_index, live_index) = parameters.doubles("live_phase_duration");
	
	Cycle_Model* apoptosis_model = cell_defaults.phenotype.death.models[apoptosis_model_index];
	int apoptotic_index = apoptosis_model->find_phase_index(PhysiCell_constants::apoptotic);
	int debris_index = apoptosis_model->find_phase_index(PhysiCell_constants::debris);
	apoptosis_model->transition_rate(apoptotic_index, debris_index) = parameters.doubles("apoptotic_duration");

	cell_defaults.phenotype.death.rates[apoptosis_model_index] = parameters.doubles("apoptotic_rate"); 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 40; 

	cell_defaults.phenotype.secretion.uptake_rates[tnf_substrate_index] = parameters.doubles("tnf_uptake_rate"); 
	cell_defaults.phenotype.secretion.secretion_rates[tnf_substrate_index] = parameters.doubles("tnf_secretion_rate");
	cell_defaults.phenotype.secretion.saturation_densities[tnf_substrate_index] = 0.5; 
	
	// add custom data here, if any 
	cell_defaults.functions.custom_cell_rule = boolean_network_rule;

	cell_defaults.phenotype.motility.migration_bias = parameters.doubles("migration_bias");
	cell_defaults.phenotype.motility.migration_speed = parameters.doubles("migration_speed");
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("cell_cell_adhesion_strength");
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("cell_cell_repulsion_strength");

	cell_defaults.parameters.max_necrosis_rate = parameters.doubles("max_necrosis_rate");
	cell_defaults.parameters.o2_necrosis_max = parameters.doubles("o2_necrosis_max");
	cell_defaults.parameters.o2_proliferation_threshold = parameters.doubles("o2_proliferation_threshold");
	//not implemented in physiboss
	//cell_defaults.parameters.o2_hypoxic_threshold = parameters.doubles("o2_hypoxic_threshold");
	cell_defaults.parameters.o2_reference = parameters.doubles("o2_reference");
	cell_defaults.parameters.o2_proliferation_saturation = parameters.doubles("o2_proliferation_saturation");
	
	cell_defaults.phenotype.geometry.radius = parameters.doubles("radius");
	cell_defaults.phenotype.geometry.nuclear_radius = parameters.doubles("nuclear_radius");
	cell_defaults.phenotype.volume.fluid_fraction = parameters.doubles("fluid_fraction");
	cell_defaults.phenotype.volume.nuclear_solid = parameters.doubles("nuclear_solid");
	cell_defaults.phenotype.volume.cytoplasmic_solid = parameters.doubles("cytoplasmic_solid");
	cell_defaults.phenotype.volume.cytoplasmic_to_nuclear_ratio = parameters.doubles("cytoplasmic_to_nuclear_ratio");

	cell_defaults.custom_data.add_variable("got_activated", "", 0);

	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/*	
	default_microenvironment_options.X_range = {-500, 500}; 
	default_microenvironment_options.Y_range = {-500, 500}; 
	default_microenvironment_options.Z_range = {-500, 500}; 
*/	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	
	
	
/*
	// all this is in XML as of August 2019 (1.6.0)
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	std::string bnd_file = parameters.strings("bnd_file");
	std::string cfg_file = parameters.strings("cfg_file");
	for (int i = 0; i < cells.size(); i++)
	{
		int x = cells[i].x;
		int y = cells[i].y;
		int z = cells[i].z;

		pC = create_cell(); 
		pC->assign_position( x, y, z );
		MaBoSSNetwork* maboss = new MaBoSSNetwork(bnd_file, cfg_file);
		pC->maboss_cycle_network = new CellCycleNetwork(maboss);
	}

	return; 
}

void boolean_network_rule(Cell* pCell, Phenotype& phenotype, double dt )
{
	if (pCell->maboss_cycle_network != NULL) {
		std::vector<bool> * nodes = pCell->maboss_cycle_network->get_nodes();
		set_input_nodes(pCell, nodes);

		pCell->maboss_cycle_network->run_maboss();

		from_nodes_to_cell(pCell, phenotype, dt);
	}
}

void set_input_nodes(Cell* pCell, std::vector<bool> * nodes) {
	
	int tnf_maboss_index = pCell->maboss_cycle_network->get_maboss()->get_node_index("TNF");
	int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 

	if (tnf_maboss_index != -1 && tnf_substrate_index != -1)
	{
		(*nodes)[tnf_maboss_index] = (*pCell->internalized_substrates)[tnf_substrate_index] > parameters.doubles("tnf_threshold");
	}
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	
	// if the cell is motile and not dead, paint it black 
	
	if( pCell->phenotype.death.dead == false && 
		pCell->type == 1 )
	{
		 output[0] = "black"; 
		 output[2] = "black"; 	
	}
	
	return output; 
}

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

		result.push_back(record);
	} while (!fin.eof());
	
	return result;
}

/* Go to proliferative if needed */
void do_proliferation( Cell* pCell, Phenotype& phenotype, double dt )
{
	//do nothing
}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
{
	MaBoSSNetwork* maboss = pCell->maboss_cycle_network->get_maboss();
	std::vector<bool>* nodes = pCell->maboss_cycle_network->get_nodes();

	int bn_index = maboss->get_node_index( "Survival" );
	if ( bn_index != -1 && (*nodes)[bn_index])
		do_proliferation( pCell, phenotype, dt );

	bn_index = maboss->get_node_index( "Apoptosis" );
	if ( bn_index != -1 && (*nodes)[bn_index] )
	{
		int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
		pCell->start_death(apoptosis_model_index);
		delete pCell->maboss_cycle_network;
		pCell->maboss_cycle_network = NULL;
		return;
	}

	bn_index = maboss->get_node_index( "NonACD" );
	if ( bn_index != -1 && (*nodes)[bn_index] )
	{
		int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
		pCell->start_death(necrosis_model_index);
		delete pCell->maboss_cycle_network;
		pCell->maboss_cycle_network = NULL;
		return;
	}

	// For model with TNF production
	bn_index = maboss->get_node_index( "NFkB" );
	if ( bn_index != -1 )
	{
		int tnf_substrate_index = microenvironment.find_density_index( "tnf" ); 
		// produce some TNF
		if ( (*nodes)[bn_index] )
		{
			pCell->phenotype.secretion.secretion_rates[tnf_substrate_index] = parameters.doubles("secretion_tnf")/microenvironment.voxels(pCell->get_current_voxel_index()).volume;
			int activated_index = pCell->custom_data.find_variable_index( "got_activated" );
			if ( pCell->custom_data[activated_index] == 0 )
				pCell->custom_data[activated_index] = 1;
		}
		else
			pCell->phenotype.secretion.secretion_rates[tnf_substrate_index] = 0;
		pCell->set_internal_uptake_constants(dt);
	}
}

void remove_density( int density_index )
{	
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ )
	{
		microenvironment.density_vector(n)[density_index] = 0; 	
	}
}

double norm(double cent[3]){ 
	return sqrt( cent[0]*cent[0] + cent[1]*cent[1] + cent[2]*cent[2] ); 
}

void inject_density(int density_index, double concentration) 
{
	// Inject given concentration on the extremities only
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ )
	{
		double cent[3] = {microenvironment.voxels(n).center[0], microenvironment.voxels(n).center[1], microenvironment.voxels(n).center[2]};

		if ( ! ((parameters.doubles("membrane_length") - norm(cent)) > 0) )
			microenvironment.density_vector(n)[density_index] = concentration; 	
	}
}