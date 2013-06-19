#ifndef CEPSCAN_HELPER_H
#define CEPSCAN_HELPER_H

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

namespace CEPscan_helper
{
	void prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	bool loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, const std::string &fileName );
	
	void streamSetParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, std::ostream &output, const std::string &prefix = "");
	
	bool checkConsistencyOfParameters( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet );
	
	void fillSetWithRange( const double min, const double max, const double step, std::set< double > &toFill );
 
	struct resultForOutput
	{
		double kappa_N;
		double lambda_N;
		double lambda_6_N;
		double yukawa_N;
		double magnetization;
		double staggered_magnetization;
		double potential;
		double d2U_ov_dmdm;
		double d2U_ov_dsds;
		double d2U_ov_dmds;
	};
	
	void findStartByTwo_1d_Scans(std::set< double > &testValues_magnetization, std::set< double > &testValues_staggeredMagnetization, constrainedEffectivePotential &CEP, double &start_m, double &start_s);
	
	void findStartByOne_2d_Scan(std::set< double > &testValues_magnetization, std::set< double > &testValues_staggeredMagnetization, constrainedEffectivePotential &CEP, double &start_m, double &start_s);
	
	bool printResultsVectorToStream(const std::vector< resultForOutput > &results, std::ostream &output);
}




#endif
