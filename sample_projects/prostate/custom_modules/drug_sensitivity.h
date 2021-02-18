#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <map>
#include <cmath>
#include <algorithm>
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

/**
 *	\main drug_AGS custom
 *	\brief Custom module file for drug_AGS example
 * 
 *	\details Modules needed for the drug_AGS example. This custom module can be used to study the inhibition of AGS cell lines with AKT, beta-catenin and TAK inhibitors.
 *
 *
 *	\date 19/10/2020
 *	\author Annika Meert, BSC-CNS
 */

using namespace std;
using namespace PhysiCell;
using namespace BioFVM;

extern const vector<pair<string, int>> cell_line_ids;
extern const vector<pair<string, int>> drug_ids;
extern const vector<pair<string, string>> drug_targets;

std::vector<std::pair<std::string, std::vector<double>>> read_csv(std::string filename);

int get_index(Cell* pCell, string drug_name, string cell_line_name);
vector<double> get_drug_sensitivity_values (Cell* pCell, string drug_name, string cell_line_name);

double get_x_from_conc(double x_conc, double max_conc);
double get_conc_from_x (double x, double max_conc);
double get_lx_from_x (double x, double max_conc);
double linear_mixed_model_function(double lx, double max_conc, double xmid, double scale);
double get_cell_viability_for_drug_conc (Cell* pCell, string cell_line, string drug_name, int index);