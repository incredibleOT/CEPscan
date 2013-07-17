// This program should simply plot the CEP. It will be able to scan in m and s and also in the other parameters.
// It will be less sopisticated (no checks for consitency etc) but will also work with inputfile

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "constrainedEffectivePotential.h"
#include "CEPscan_helper.h"
#include "plotPotential_helper.h"

using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
	std::map< std::string, double > parametersDouble;
	std::map< std::string, int > parametersInt;
	std::map< std::string, std::string > parametersString;
	std::map< std::string, bool > parametersIsSet;
	
	plotPotential_helper::prepareParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet);
	size_t numberOfParameters=parametersDouble.size()+parametersInt.size()+parametersString.size();
	
	if(narg!=2)
	{
		cerr <<"Error, start program with:" <<endl;
		cerr <<arg[0] <<"   inputFile" <<endl;
		exit(EXIT_FAILURE);
	}
	
	if( !CEPscan_helper::loadParameterMapsFromFile(parametersDouble, parametersInt, parametersString, parametersIsSet, arg[1]) )
	{
		cerr <<"Error loading input file" <<endl <<arg[1] <<endl;
		exit(EXIT_FAILURE);
	}
	
	cout <<endl <<"Parameters loaded:" <<endl;
	CEPscan_helper::streamSetParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet,cout);
	
	int Ls(parametersInt["L_s"]), Lt(parametersInt["L_t"]);
	constrainedEffectivePotential CEP(Ls,Ls,Ls,Lt,parametersInt["antiperiodic_L_t"]);
	
	//set N_f, rho, r if set an different
	if(parametersIsSet["N_f"] && parametersInt["N_f"]!=CEP.get_N_f()){ CEP.set_N_f(parametersInt["N_f"]); }
	if(parametersIsSet["rho"] && parametersDouble["rho"]!=CEP.get_rho()){ CEP.set_rho(parametersDouble["rho"]); }
	if(parametersIsSet["r"] && parametersDouble["r"]!=CEP.get_r()){ CEP.set_rho(parametersDouble["r"]); }
	
	//include_bosonic_loop
	CEP.set_useBosonicLoop(parametersInt["include_bosonic_loop"]);
	//use_improved_treeLevel
	CEP.set_useImprovedGaussian(parametersInt["use_improved_treeLevel"]);
	//use_improved_firstOrder
	CEP.set_useImprovedFirstOrder(parametersInt["use_improved_firstOrder"]);
	
	std::set< double > kappa_N_values, lambda_N_values, lambda_6_N_values, yukawa_N_values, mag_values, stag_values; //for parameters
	std::set< double > testValues_magnetization, testValues_staggeredMagnetization; //for scanning of startingpoint
	
	//kappa_N
	if(parametersIsSet["scan_kappa_N"] && parametersInt["scan_kappa_N"] )
	{
		CEPscan_helper::fillSetWithRange( parametersDouble["kappa_N_min"], parametersDouble["kappa_N_max"], parametersDouble["kappa_N_step"], kappa_N_values );
	}
	else{ kappa_N_values.insert(parametersDouble["kappa_N"]); }
	//lambda_N
	if(parametersIsSet["scan_lambda_N"] && parametersInt["scan_lambda_N"] )
	{
		CEPscan_helper::fillSetWithRange( parametersDouble["lambda_N_min"], parametersDouble["lambda_N_max"], parametersDouble["lambda_N_step"], lambda_N_values );
	}
	else{ lambda_N_values.insert(parametersDouble["lambda_N"]); }
	//lambda_6_N
	if(parametersIsSet["scan_lambda_6_N"] && parametersInt["scan_lambda_6_N"] )
	{
		CEPscan_helper::fillSetWithRange( parametersDouble["lambda_6_N_min"], parametersDouble["lambda_6_N_max"], parametersDouble["lambda_6_N_step"], lambda_6_N_values );
	}
	else{ lambda_6_N_values.insert(parametersDouble["lambda_6_N"]); }
	//yukawa_N
	if(parametersIsSet["scan_yukawa_N"] && parametersInt["scan_yukawa_N"] )
	{
		CEPscan_helper::fillSetWithRange( parametersDouble["yukawa_N_min"], parametersDouble["yukawa_N_max"], parametersDouble["yukawa_N_step"], yukawa_N_values );
	}
	else{ yukawa_N_values.insert(parametersDouble["yukawa_N"]); }
	//mag
	if(parametersIsSet["scan_mag"] && parametersInt["scan_mag"] )
	{
		CEPscan_helper::fillSetWithRange( parametersDouble["mag_min"], parametersDouble["mag_max"], parametersDouble["mag_step"], mag_values );
	}
	else{ mag_values.insert(parametersDouble["mag"]); }
	//stag
	if(parametersIsSet["scan_stag"] && parametersInt["scan_stag"] )
	{
		CEPscan_helper::fillSetWithRange( parametersDouble["stag_min"], parametersDouble["stag_max"], parametersDouble["stag_step"], stag_values );
	}
	else{ stag_values.insert(parametersDouble["stag"]); }
	
	typedef plotPotential_helper::resultForOutput resultType;
	std::vector< resultType > results;
	cout <<"start scanning" <<endl;
	//now iterate
	for(std::set< double >::const_iterator kappa_N=kappa_N_values.begin(); kappa_N!=kappa_N_values.end(); ++kappa_N)
	{
		cout <<"k=" <<*kappa_N <<endl;
		double lambdaFactor=(parametersInt["interpret_lambda_as_continuum"])?( 4.0 * (*kappa_N) * (*kappa_N) ) : 1.0;
		double lambda_6_Factor=(parametersInt["interpret_lambda_6_as_continuum"])?( 8.0 * (*kappa_N) * (*kappa_N) * (*kappa_N) ) : 1.0;
		double yukawaFactor=(parametersInt["interpret_yukawa_as_continuum"])?( sqrt(2.0* (*kappa_N)) ): 1.0;
		for(std::set< double >::const_iterator lambda_N=lambda_N_values.begin(); lambda_N!=lambda_N_values.end(); ++lambda_N)
		{
			cout <<"    l=" <<*lambda_N <<endl;
			for(std::set< double >::const_iterator lambda_6_N=lambda_6_N_values.begin(); lambda_6_N!=lambda_6_N_values.end(); ++lambda_6_N)
			{
				cout <<"        l_6="<<*lambda_6_N <<endl;
				for(std::set< double >::const_iterator yukawa_N=yukawa_N_values.begin(); yukawa_N!=yukawa_N_values.end(); ++yukawa_N)
				{
					cout <<"            y=" <<*yukawa_N <<endl;
					CEP.set_kappa_lambda_lambda_6_yukawa_N( *kappa_N, *lambda_N * lambdaFactor, *lambda_6_N * lambda_6_Factor, *yukawa_N * yukawaFactor);
					for(std::set< double >::const_iterator mag=mag_values.begin(); mag!=mag_values.end(); ++mag)
					{
						cout <<"                m=" <<*mag;
						for(std::set< double >::const_iterator stag=stag_values.begin(); stag!=stag_values.end(); ++stag)
						{
							cout <<"    s=" <<*stag <<endl;
							resultType dummyResult;
							
							CEP.computeConstrainedEffectivePotential_FunctionAndGradient(*mag, *stag, dummyResult.potential, dummyResult.dU_ov_dm, dummyResult.dU_ov_ds);
							CEP.computeConstrainedEffectivePotential_secondDerivatives(*mag, *stag, dummyResult.d2U_ov_dmdm, dummyResult.d2U_ov_dsds, dummyResult.d2U_ov_dmds);
							
							dummyResult.kappa_N=CEP.get_kappa_N();
							dummyResult.lambda_N=CEP.get_lambda_N();
							dummyResult.lambda_6_N=CEP.get_lambda_6_N();
							dummyResult.yukawa_N=CEP.get_yukawa_N();
							dummyResult.magnetization=*mag;
							dummyResult.staggered_magnetization=*stag;
							
							results.push_back( dummyResult );
						}//stag
					}//mag
				}//yukawa_N
			}//lambda_6_N
		}//lambda_N
	}//kappa_N
	
	{
		bool fileOK(false);
		std::ofstream outputFile;
		if(parametersIsSet["outputfile"])
		{
			std::string orig( parametersString["outputfile"] );
			std::string outputFileName( plotPotential_helper::generateOutputFileName( orig, parametersDouble, parametersInt, parametersIsSet) );
			
			cout <<"output fileName: " << outputFileName <<endl;
			
			outputFile.open( outputFileName.c_str() );
			if(!outputFile.good())
			{
				cerr <<"Error opening output file" <<endl <<outputFileName <<endl;
				cerr <<"send ouput to cout!" <<endl;
			}
			else
			{
				cout <<"opened file:" <<endl <<outputFileName <<endl;
				fileOK=true; 
			}
		}
		std::ostream &output=fileOK?outputFile:cout;
		
		output <<"# Output of CEPscan - Plotting the constrained effective Potential" <<endl;
		output <<"# parameters set:" <<endl;
		CEPscan_helper::streamSetParameterMaps( parametersDouble, parametersInt, parametersString, parametersIsSet, outputFile, "#");
		output <<"# Output format is:" <<endl;
		output <<"# kappa_N   lambda_N  lambda_6_N  yukawa_N   Magnetization   stageredMagnetization   U_min   dU_ov_dm   dU_ov_ds   d2U_ov_dmdm   d2U_ov_dsds   d2U_ov_dmds" <<endl;
		plotPotential_helper::printResultsVectorToStream( results, output );
		outputFile.close();
	}
			
			
			
			
			
			
	
	if(parametersDouble.size()+parametersInt.size()+parametersString.size() != numberOfParameters || parametersIsSet.size()!=numberOfParameters )
	{
		cerr <<"Error, number of parameters changed!" <<endl;
		cerr <<"parametersDouble.size()=" <<parametersDouble.size() <<endl;
		cerr <<"parametersInt.size()=" <<parametersInt.size() <<endl;
		cerr <<"parametersString.size()=" <<parametersString.size() <<endl;
		cerr <<"parametersIsSet.size()=" <<parametersIsSet.size() <<endl;
		exit(EXIT_FAILURE);
	}
}
	
	