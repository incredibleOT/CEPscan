#ifndef PLOTPOTENTIAL_HELPER_H
#define PLOTPOTENTIAL_HELPER_H

#include <istream>
#include <iostream>
#include <fstream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "constrainedEffectivePotential.h"

namespace plotPotential_helper
{
	void prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	struct resultForOutput
	{
		double kappa_N;
		double lambda_N;
		double lambda_6_N;
		double yukawa_N;
		double magnetization;
		double staggered_magnetization;
		double potential;
		double dU_ov_dm;
		double dU_ov_ds;
		double d2U_ov_dmdm;
		double d2U_ov_dsds;
		double d2U_ov_dmds;
	};
	
	std::string generateOutputFileName( const std::string &original, std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, bool > &paraIsSet );
	
	bool printResultsVectorToStream(const std::vector< resultForOutput > &results, std::ostream &output);
}	
	
	
#endif