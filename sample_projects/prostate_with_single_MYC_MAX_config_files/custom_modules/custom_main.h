
#include "../core/PhysiCell.h"

/**
 *	\main drug_AGS custom main file
 *	\brief Custom module file for drug_AGS example
 * 
 *	\details Modules needed for the drug_AGS example. This custom module can be used to study the inhibition of AGS cell lines with AKT, beta-catenin and TAK inhibitors.
 *
 *
 *	\date 19/10/2020
 *	\author Arnau Montagud, BSC-CNS, with code previously developed by Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

using namespace BioFVM; 

void inject_density_sphere(int density_index, double concentration, double membrane_lenght);
void remove_density( int density_index );
void set_density_for_current_time (int density_index, double current_time, double max_time, double time_add_dens, double time_put_dens, double duration_add_dens, double time_remove_dens, double time_dens_next, double concentration_dens, double membrane_length );
inline double current_value( double min, double max, double percent )
{ return (min + (max-min) * percent); };

static const double EPSILON = std::numeric_limits<double>::epsilon();

/** \brief Relative difference between two numbers */
inline double relative_diff( double a, double b )
{ if ( b < EPSILON ) return 0; return ( fabs(a-b)/b ); };


/* Change the current value of the input coefficient, increase or decrease according to up value */
void evolve_coef( int up, double* coef, double dt );