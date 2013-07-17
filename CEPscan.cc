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
	
	//include_bosonic_loop
	CEP.set_useBosonicLoop(parametersInt["include_bosonic_loop"]);
	//use_improved_treeLevel
	CEP.set_useImprovedGaussian(parametersInt["use_improved_treeLevel"]);
	//use_improved_firstOrder
	CEP.set_useImprovedFirstOrder(parametersInt["use_improved_firstOrder"]);
	
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
	std::set< double > kappa_N_values, lambda_N_values, lambda_6_N_values, yukawa_N_values; //for parameters
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
	
	
	//set starting values
	CEP.set_kappa_N(*kappa_N_values.begin());
	CEP.set_lambda_N(*lambda_N_values.begin());
	CEP.set_lambda_6_N(*lambda_6_N_values.begin());
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
	cout <<"start scanning" <<endl;
	//now iterate
	for(std::set< double >::const_iterator kappa_N=kappa_N_values.begin(); kappa_N!=kappa_N_values.end(); ++kappa_N)
	{
// 		cout <<"kappa=" <<*kappa_N <<endl;
		double lambdaFactor=(parametersInt["interpret_lambda_as_continuum"])?( 4.0 * (*kappa_N) * (*kappa_N) ) : 1.0;
		double lambda_6_Factor=(parametersInt["interpret_lambda_6_as_continuum"])?( 8.0 * (*kappa_N) * (*kappa_N) * (*kappa_N) ) : 1.0;
		double yukawaFactor=(parametersInt["interpret_yukawa_as_continuum"])?( sqrt(2.0* (*kappa_N)) ): 1.0;
		for(std::set< double >::const_iterator lambda_N=lambda_N_values.begin(); lambda_N!=lambda_N_values.end(); ++lambda_N)
		{
// 			cout <<"lambda=" <<*lambda_N <<endl;
			for(std::set< double >::const_iterator lambda_6_N=lambda_6_N_values.begin(); lambda_6_N!=lambda_6_N_values.end(); ++lambda_6_N)
			{
				for(std::set< double >::const_iterator yukawa_N=yukawa_N_values.begin(); yukawa_N!=yukawa_N_values.end(); ++yukawa_N)
				{
					CEP.set_kappa_lambda_lambda_6_yukawa_N( *kappa_N, *lambda_N * lambdaFactor, *lambda_6_N * lambda_6_Factor, *yukawa_N * yukawaFactor);
						
					cout <<"kappa_N=" <<CEP.get_kappa_N() <<"  lambda_N=" <<CEP.get_lambda_N() <<"  lambda_6_N=" <<CEP.get_lambda_6_N() <<"  yukawa_N=" <<CEP.get_yukawa_N(); //<<endl
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
						dummyResult.lambda_6_N=CEP.get_lambda_6_N();
						dummyResult.yukawa_N=CEP.get_yukawa_N();
						dummyResult.magnetization=std::abs(start_m);
						dummyResult.staggered_magnetization=std::abs(start_s);
						dummyResult.potential=CEP.getActualMinimizerValue();
						CEP.getActualSecondDerivative(dummyResult.d2U_ov_dmdm, dummyResult.d2U_ov_dsds, dummyResult.d2U_ov_dmds);
						results.push_back(dummyResult);
					}
						
				}//yukawa
			}//lambda_6
		}//lambda
	}//kappa
	
	
	//output to cout
	{
		cout <<"Results (kappa_N,   lambda_N,   lambda_6_N,   yukawa_N,   mag,   stag.mag,   pot,   d2U_ov_dmdm,    d2U_ov_dsds,    d2U_ov_dmds:" <<endl;
		CEPscan_helper::printResultsVectorToStream( results, cout );
	}
	//
	{
		bool fileOK(false);
		std::ofstream outputFile;
		if(parametersIsSet["outputfile"])
		{
			std::string outputFileName(parametersString["outputfile"]);
			if( outputFileName.find("[Ls]")!=std::string::npos )
			{
				std::ostringstream ss;
				ss <<"L" <<parametersInt["L_s"];
				outputFileName.replace(outputFileName.find("[Ls]"),4, ss.str() );
			}
			if( outputFileName.find("[Lt]")!=std::string::npos )
			{
				std::ostringstream ss;
				ss <<"T" <<parametersInt["L_t"];
				outputFileName.replace(outputFileName.find("[Lt]"),4, ss.str() );
			}
			if( outputFileName.find("[k]")!=std::string::npos )
			{
				std::ostringstream ss;
				ss <<"k_";
				ss.precision(7);
				if(parametersInt["scan_kappa_N"]){ss <<parametersDouble["kappa_N_min"] <<"_"<<parametersDouble["kappa_N_max"]; }
				else {ss <<parametersDouble["kappa_N"]; }
				outputFileName.replace(outputFileName.find("[k]"),3, ss.str() );
			}
			if( outputFileName.find("[l]")!=std::string::npos )
			{
				std::ostringstream ss;
				ss <<"l_";
				ss.precision(7);
				if(parametersInt["scan_lambda_N"]){ss <<parametersDouble["lambda_N_min"] <<"_"<<parametersDouble["lambda_N_max"]; }
				else {ss <<parametersDouble["lambda_N"]; }
				outputFileName.replace(outputFileName.find("[l]"),3, ss.str() );
			}
			if( outputFileName.find("[l6]")!=std::string::npos )
			{
				std::ostringstream ss;
				ss <<"l6_";
				ss.precision(7);
				if(parametersInt["scan_lambda_6_N"]){ss <<parametersDouble["lambda_6_N_min"] <<"_"<<parametersDouble["lambda_6_N_max"]; }
				else {ss <<parametersDouble["lambda_6_N"]; }
				outputFileName.replace(outputFileName.find("[l6]"),4, ss.str() );
			}
			if( outputFileName.find("[y]")!=std::string::npos )
			{
				std::ostringstream ss;
				ss <<"y_";
				ss.precision(7);
				if(parametersInt["scan_yukawa_N"]){ss <<parametersDouble["yukawa_N_min"] <<"_"<<parametersDouble["yukawa_N_max"]; }
				else {ss <<parametersDouble["yukawa_N"]; }
				outputFileName.replace(outputFileName.find("[y]"),3, ss.str() );
			}
			if( outputFileName.find("[loop]")!=std::string::npos )
			{
				std::ostringstream ss;
				if(parametersInt["include_bosonic_loop"]){ ss <<"withLoop";}
				else{ ss <<"noLoop";}
				outputFileName.replace(outputFileName.find("[loop]"),6, ss.str() );
			}
			if( outputFileName.find("[bc]")!=std::string::npos )
			{
				std::ostringstream ss;
				if(parametersInt["antiperiodic_L_t"]){ ss <<"aBC";}
				else{ ss <<"pBC";}
				outputFileName.replace(outputFileName.find("[bc]"),4, ss.str() );
			}
			if( outputFileName.find("[det]")!=std::string::npos )
			{
				std::ostringstream ss;
				if(parametersInt["use_improved_treeLevel"]){ ss <<"withDet";}
				else{ ss <<"noDet";}
				outputFileName.replace(outputFileName.find("[det]"),5, ss.str() );
			}
			if( outputFileName.find("[1st]")!=std::string::npos )
			{
				std::ostringstream ss;
				if(parametersInt["use_improved_firstOrder"]){ ss <<"with1st";}
				else{ ss <<"no1st";}
				outputFileName.replace(outputFileName.find("[1st]"),5, ss.str() );
			}
			
			
			
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
		
		output <<"# Output of CEPscan - Minimization of the constrained effective Potential" <<endl;
		output <<"# parameters set:" <<endl;
		CEPscan_helper::streamSetParameterMaps( parametersDouble, parametersInt, parametersString, parametersIsSet, outputFile, "#");
		output <<"# Output format is:" <<endl;
		output <<"# kappa_N   lambda_N  lambda_6_N  yukawa_N   Magnetization   stageredMagnetization   U_min   d2U_ov_dmdm   d2U_ov_dsds   d2U_ov_dmds" <<endl;
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

