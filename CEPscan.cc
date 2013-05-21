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

using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
	std::map< std::string, double > parametersDouble;
	std::map< std::string, int > parametersInt;
	std::map< std::string, std::string > parametersString;
	std::map< std::string, bool > parametersIsSet;
	
	CEPscan_helper::prepareParameterMaps(parametersDouble, parametersInt, parametersString, parametersIsSet);
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
	
	cout <<endl <<"Check consistency of parameters" <<endl;
	if(!(CEPscan_helper::checkConsistencyOfParameters(parametersDouble, parametersInt, parametersString, parametersIsSet)))
	{
		cerr <<"Failed! Some parameters are inconsistent!" <<endl;
		exit(EXIT_FAILURE);
	}
	cout <<"passed" <<endl;
	
	
	
	//initiallize CEP
	int Ls(parametersInt["L_s"]), Lt(parametersInt["L_t"]);
	constrainedEffectivePotential CEP(Ls,Ls,Ls,Lt,parametersInt["antiperiodic_L_t"]);
	
	//set N_f, rho, r if set an different
	if(parametersIsSet["N_f"] && parametersInt["N_f"]!=CEP.get_N_f()){ CEP.set_N_f(parametersInt["N_f"]); }
	if(parametersIsSet["rho"] && parametersDouble["rho"]!=CEP.get_rho()){ CEP.set_rho(parametersDouble["rho"]); }
	if(parametersIsSet["r"] && parametersDouble["r"]!=CEP.get_r()){ CEP.set_rho(parametersDouble["r"]); }
	
	//set minimization details
	CEP.set_toleranceForLineMinimization(parametersDouble["tolerance_for_line_minimization"]);
	if(parametersIsSet["tolerance_for_convergence"])
	{
		CEP.set_toleranceForConvergence(parametersDouble["tolerance_for_convergence"]);
	}
	switch(parametersInt["starting_procedure"])
	{
		case 1: //same as 2
		case 2: 
		CEP.set_initialStepSize(parametersDouble["initial_step_size"]);
		break;
		case 3://same as 4
		case 4:
		CEP.set_initialStepSize(parametersDouble["step_m"]<=parametersDouble["step_s"] ? 0.05*parametersDouble["step_m"] : 0.05*parametersDouble["step_s"]);
		break;
		default: 
		cerr <<"Error, no starting_procedure=" <<parametersInt["starting_procedure"] <<endl;
		exit(EXIT_FAILURE);
	}
	CEP.set_maxNumerOfIterations(parametersInt["max_numer_of_iterations"]);
	CEP.set_minimizationAlgorithm(parametersInt["minimization_algorithm"]);
	
	//now prepare the scanning lists
	std::set< double > kappa_N_values, lambda_N_values, yukawa_N_values; //for aprameters
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
	//yukawa_N
	if(parametersIsSet["scan_yukawa_N"] && parametersInt["scan_yukawa_N"] )
	{
		CEPscan_helper::fillSetWithRange( parametersDouble["yukawa_N_min"], parametersDouble["yukawa_N_max"], parametersDouble["yukawa_N_step"], yukawa_N_values );
	}
	else{ yukawa_N_values.insert(parametersDouble["yukawa_N"]); }
	
	
	//set starting values
	CEP.set_kappa_N(*kappa_N_values.begin());
	CEP.set_lambda_N(*lambda_N_values.begin());
	CEP.set_yukawa_N(*yukawa_N_values.begin());
	
	if(parametersInt["starting_procedure"]==1 || parametersInt["starting_procedure"]==2)
	{
		CEP.initializeMinimizer(parametersDouble["start_m"] , parametersDouble["start_s"]);
	}
	else
	{
		CEP.initializeMinimizer(1.0, 1.0);
		//prepare startingValues
		CEPscan_helper::fillSetWithRange( parametersDouble["minimal_m"], parametersDouble["maximal_m"], parametersDouble["step_m"], testValues_magnetization);
		CEPscan_helper::fillSetWithRange( parametersDouble["minimal_s"], parametersDouble["maximal_s"], parametersDouble["step_s"], testValues_staggeredMagnetization);
	}
	
	typedef CEPscan_helper::resultForOutput resultType;
	std::vector< resultType > results;
	
	//now iterate
	for(std::set< double >::const_iterator kappa_N=kappa_N_values.begin(); kappa_N!=kappa_N_values.end(); ++kappa_N)
	{
		double lambdaFactor=(parametersInt["interpret_lambda_as_continuum"])?( 4.0 * (*kappa_N) * (*kappa_N) ) : 1.0;
		double yukawaFactor=(parametersInt["interpret_yukawa_as_continuum"])?( sqrt(2.0* (*kappa_N)) ): 1.0;
		for(std::set< double >::const_iterator lambda_N=lambda_N_values.begin(); lambda_N!=lambda_N_values.end(); ++lambda_N)
		{
			for(std::set< double >::const_iterator yukawa_N=yukawa_N_values.begin(); yukawa_N!=yukawa_N_values.end(); ++yukawa_N)
			{
				CEP.set_kappa_lambda_yukawa_N( *kappa_N, *lambda_N * lambdaFactor, *yukawa_N * yukawaFactor);
					
				cout <<"kappa_N=" <<CEP.get_kappa_N() <<"  lambda_N=" <<CEP.get_lambda_N() <<"  yukawa_N=" <<CEP.get_yukawa_N(); //<<endl
				double start_m(0.0), start_s(0.0);
				int nIter;
				switch( parametersInt["starting_procedure"] )
				{
					case 1:
					CEP.reInitializeMinimizer(parametersDouble["start_m"] , parametersDouble["start_s"]); 
					break;
					case 2: break;
					case 3: 
					CEPscan_helper::findStartByTwo_1d_Scans(testValues_magnetization, testValues_staggeredMagnetization, CEP, start_m, start_s);
					CEP.reInitializeMinimizer(start_m, start_s);
					break;
					case 4:
					CEPscan_helper::findStartByOne_2d_Scan(testValues_magnetization, testValues_staggeredMagnetization, CEP, start_m, start_s);
					CEP.reInitializeMinimizer(start_m, start_s);
					break;
					default: 
					cerr <<"Error, no starting_procedure=" <<parametersInt["starting_procedure"] <<endl;
					exit(EXIT_FAILURE);
				}
				CEP.getActualMinimizerLocation(start_m, start_s);
				cout <<"  start at:  m=" <<start_m <<"  s=" <<start_s;// <<endl;
				switch( parametersInt["iteration_stopping_criteria"] )
				{
					case 1: 
					nIter=CEP.itarateUntilToleranceReached(parametersDouble["tolerance_for_convergence"]);
					break;
					case 2:
					nIter=CEP.iterateUntilIterationStopps();
					break;
					default:
					cerr <<"Error, no iteration_stopping_criteria=" <<parametersInt["iteration_stopping_criteria"] <<endl;
					exit(EXIT_FAILURE);
				}
				if(nIter==0)
				{
					cout <<endl <<"Minimization did not converge in " <<CEP.get_maxNumerOfIterations() <<" steps" <<endl;
				}
				else
				{
					CEP.getActualMinimizerLocation(start_m, start_s);
					cout <<"  Minimum: U(m=" <<start_m <<"  s=" <<start_s <<")=" <<CEP.getActualMinimizerValue();
					cout <<"  (" <<nIter <<" iter.)" <<endl;
					resultType dummyResult;
					dummyResult.kappa_N=CEP.get_kappa_N();
					dummyResult.lambda_N=CEP.get_lambda_N();
					dummyResult.yukawa_N=CEP.get_yukawa_N();
					dummyResult.magnetization=std::abs(start_m);
					dummyResult.staggered_magnetization=std::abs(start_s);
					dummyResult.potential=CEP.getActualMinimizerValue();
					results.push_back(dummyResult);
				}
					
			}//yukawa
		}//lambda
	}//kappa
	
	
	//output to cout
	{
		cout <<"Results (kappa_N,   lambda_N,   yukawa_N,   mag,   stag.mag,   pot:" <<endl;
		CEPscan_helper::printResultsVectorToStream( results, cout );
	}
	//
	{
		bool fileOK(false);
		std::ofstream outputFile;
		if(parametersIsSet["outputfile"])
		{
			outputFile.open( parametersString["outputfile"].c_str() );
			if(!outputFile.good())
			{
				cerr <<"Error opening output file" <<endl <<parametersString["outputfile"] <<endl;
				cerr <<"send ouput to cout!" <<endl;
			}
			else
			{
				cout <<"opened file:" <<endl <<parametersString["outputfile"] <<endl;
				fileOK=true; 
			}
		}
		std::ostream &output=fileOK?outputFile:cout;
		
		output <<"# Output of CEPscan - Minimization of the constrained effective Potential" <<endl;
		output <<"# parameters set:" <<endl;
		CEPscan_helper::streamSetParameterMaps( parametersDouble, parametersInt, parametersString, parametersIsSet, outputFile, "#");
		output <<"# Output format is:" <<endl;
		output <<"# kappa_N    lambda_N    yukawa_N    Magnetization    stageredMagnetization    minimalValue" <<endl;
		CEPscan_helper::printResultsVectorToStream( results, output );
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

