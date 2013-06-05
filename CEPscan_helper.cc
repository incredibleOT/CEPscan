#include "CEPscan_helper.h"


void CEPscan_helper::prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet )
{
	paraD.clear(); paraI.clear(); paraS.clear(); paraIsSet.clear();
	
	paraI["L_s"]=-1;
	paraI["L_t"]=-1;
	paraI["antiperiodic_L_t"]=-1;
	
	paraI["scan_kappa_N"]=0;
	paraD["kappa_N"] = 0.0;
	paraD["kappa_N_min"] = 0.0;
	paraD["kappa_N_max"] = 0.0;
	paraD["kappa_N_step"] = 0.0;

	paraI["scan_lambda_N"]=0;
	paraD["lambda_N"] = 0.0;
	paraD["lambda_N_min"] = 0.0;
	paraD["lambda_N_max"] = 0.0;
	paraD["lambda_N_step"] = 0.0;
	
	paraI["scan_yukawa_N"]=0;
	paraD["yukawa_N"] = 0.0;
	paraD["yukawa_N_min"] = 0.0;
	paraD["yukawa_N_max"] = 0.0;
	paraD["yukawa_N_step"] = 0.0;
	
	paraI["interpret_yukawa_as_continuum"]=0;
	paraI["interpret_lambda_as_continuum"]=0;
	
	paraI["include_bosonic_loop"]=0;
	
	//directly set to default values
	paraI["N_f"]=1;
	paraD["rho"]=1.0;
	paraD["r"]=0.5;
	
	//1 when tolerance_for_convergence is reached,   2 when iteration stopped, since it cannot improve)
	paraI["iteration_stopping_criteria"]=-1;
	
	paraD["tolerance_for_line_minimization"]=-1.0;
	paraD["tolerance_for_convergence"]=-1.0;
	paraI["max_numer_of_iterations"]=100;
	
	
	//1 conjugate_fr, 2 conjugate_pr, 3 vector_bfgs2, 4 vector_bfgs, 5 steepest_descent
	paraI["minimization_algorithm"]=-1;
	//1 use given starting point, 2 use given starting point first and during scan last result
	//3 perform a scan in {m_0:m_1:m_s , s_0} and {m_0, s_0:s_1:s_s}, determin the lowest point and start from there - initial stepsize will be the minimum of (s_0/2,m_0/2)
	//4perform a 2-d scan in {m_0:m_1:m_s ,s_0:s_1:s_s}, determin the lowest point and start from there - initial stepsize will be the minimum of (s_0/2,m_0/2))
	paraI["starting_procedure"]=-1;
	paraD["initial_step_size"]=-1.0;
	
	paraD["start_m"]=-1.0;
	paraD["minimal_m"]=-1.0; 
	paraD["maximal_m"]=-1.0; 
	paraD["step_m"]=-1.0;    

	paraD["start_s"]=-1.0;
	paraD["minimal_s"]=-1.0;
	paraD["maximal_s"]=-1.0;
	paraD["step_s"]=-1.0;
 
	paraS["outputfile"]="";
	
	for( std::map< std::string, double >::const_iterator iter=paraD.begin(); iter!=paraD.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
	for( std::map< std::string, int >::const_iterator iter=paraI.begin(); iter!=paraI.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
	for( std::map< std::string, std::string >::const_iterator iter=paraS.begin(); iter!=paraS.end(); ++iter )
	{
		paraIsSet[iter->first]=false;
	}
}//prepareParameterMaps





bool CEPscan_helper::loadParameterMapsFromFile( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, const std::string &fileName )
{
	std::ifstream inputFile(fileName.c_str());
	std::string line,word;
	if(!(inputFile.good()))
	{
		std::cerr <<"Error opening inputfile " <<fileName <<std::endl;
		return false;
	}
	while(inputFile)
	{
		getline(inputFile, line);
		if(line.size()==0 || line[0] =='#' || line.find_first_not_of(' ') == std::string::npos){ continue; }
		std::istringstream strm(line);
		if(!(strm >> word)){ continue; }
		if( paraI.find(word) != paraI.end() )
		{
			if(paraIsSet[word]){ std::cerr<<"Error, parameter " <<word <<" already set" <<std::endl; return false; } 
			if(!(strm >> paraI[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else if( paraD.find(word) != paraD.end() )
		{ 
			if(paraIsSet[word]){ std::cerr<<"Error, parameter " <<word <<" already set" <<std::endl; return false; }
			if(!(strm >> paraD[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else if( paraS.find(word) != paraS.end() )
		{ 
			if(paraIsSet[word]){ std::cerr<<"Error, parameter " <<word <<" already set" <<std::endl; return false; }
			if(!(strm >> paraS[word]))
			{
				std::cerr <<"Error loading value of " <<word <<std::endl; 
				std::cerr <<"actual line: " <<std::endl << line <<std::endl; return false; 
			}
			else{ paraIsSet[word]=true;}
		}
		else{ std::cerr <<"Warning, there is no parameter called \"" <<word <<"\""<<std::endl; return false;}
		if(inputFile.eof()){ break; }
	}
	inputFile.close(); inputFile.clear();
	return true;
}//loadParameterMapsFromFile





void CEPscan_helper::streamSetParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet, std::ostream &output, const std::string &prefix)
{
	for( std::map< std::string, double >::const_iterator iter=paraD.begin(); iter!=paraD.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
	for( std::map< std::string, int >::const_iterator iter=paraI.begin(); iter!=paraI.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
	for( std::map< std::string, std::string >::const_iterator iter=paraS.begin(); iter!=paraS.end(); ++iter )
	{
		if(paraIsSet[iter->first])
		{
			if(!(prefix=="")){output <<prefix <<" ";}  
			output <<iter->first <<"     " <<iter->second <<std::endl;
		}
	}
}//streamSetParameterMaps




bool CEPscan_helper::checkConsistencyOfParameters( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet )
{
	//extend related
	if(!(paraIsSet["L_s"] && paraIsSet["L_t"]) || paraI["L_s"]<=0 || paraI["L_t"]<=0)
	{
		std::cerr <<"Error, no or non-positive lattice extends given" <<std::endl;
		return false;
	}
	if(paraI["L_s"]%2!=0 || paraI["L_t"]%2!=0)
	{
		std::cerr <<"Error, only even lattice extend allowed" <<std::endl;
		return false;
	}
	if(!paraIsSet["antiperiodic_L_t"])
	{
		std::cerr <<"Error, antiperiodic_L_t is not specified" <<std::endl;
		return false;
	}
	
	//kappa_N
	if( (!paraIsSet["scan_kappa_N"] || paraI["scan_kappa_N"]==0) && !paraIsSet["kappa_N"])
	{
		std::cerr <<"Error, no kappa_N given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_kappa_N"] && paraI["scan_kappa_N"] )
	{
		if( !( paraIsSet["kappa_N_min"] && paraIsSet["kappa_N_max"] && paraIsSet["kappa_N_step"] ) )
		{
			std::cerr <<"Error, no scan range in kappa_N given" <<std::endl;
			return false;
		}
		if( paraD["kappa_N_max"] < paraD["kappa_N_min"] || paraD["kappa_N_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in kappa_N" <<std::endl;
			return false;
		}
	}
	
	//lambda_N
	if( (!paraIsSet["scan_lambda_N"] || paraI["scan_lambda_N"]==0) && !paraIsSet["lambda_N"])
	{
		std::cerr <<"Error, no lambda_N given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_lambda_N"] && paraI["scan_lambda_N"] )
	{
		if( !( paraIsSet["lambda_N_min"] && paraIsSet["lambda_N_max"] && paraIsSet["lambda_N_step"] ) )
		{
			std::cerr <<"Error, no scan range in lambda_N given" <<std::endl;
			return false;
		}
		if( paraD["lambda_N_max"] < paraD["lambda_N_min"] || paraD["lambda_N_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in lambda_N" <<std::endl;
			return false;
		}
	}
	
	//yukawa_N
	if( (!paraIsSet["scan_yukawa_N"] || paraI["scan_yukawa_N"]==0) && !paraIsSet["yukawa_N"])
	{
		std::cerr <<"Error, no yukawa_N given" <<std::endl;
		return false;
	}
	if( paraIsSet["scan_yukawa_N"] && paraI["scan_yukawa_N"] )
	{
		if( !( paraIsSet["yukawa_N_min"] && paraIsSet["yukawa_N_max"] && paraIsSet["yukawa_N_step"] ) )
		{
			std::cerr <<"Error, no scan range in yukawa_N given" <<std::endl;
			return false;
		}
		if( paraD["yukawa_N_max"] < paraD["yukawa_N_min"] || paraD["yukawa_N_step"] <=0.0 )
		{
			std::cerr <<"Error, inconsistent scan range in yukawa_N" <<std::endl;
			return false;
		}
	}
	
	//interprete parameters as physical
	if( paraIsSet["interpret_yukawa_as_continuum"] && paraI["interpret_yukawa_as_continuum"] )
	{
		if( ( paraIsSet["scan_kappa_N"] && paraI["scan_kappa_N"] && paraD["kappa_N_min"]<0.0 ) || 
		    ( ( !paraIsSet["scan_kappa_N"] || !paraI["scan_kappa_N"] ) && paraD["kappa_N"]<0.0 ) )
		{
			std::cerr <<"Error, interpret_yukawa_as_continuum incompatible with negative kappa_N" <<std::endl;
			return false;
		}
	}
	
	//include_bosonic_loop
	if( !paraIsSet["include_bosonic_loop"] )
	{
		std::cerr <<"Error, include_bosonic_loop not set" <<std::endl;
		return false;
	}
	
	//iteration stopping
	if( !paraIsSet["iteration_stopping_criteria"] )
	{
		std::cerr <<"Error, no iteration_stopping_criteria set" <<std::endl;
		return false;
	}
	switch(paraI["iteration_stopping_criteria"])
	{
		case 1:
		if(!paraIsSet["tolerance_for_convergence"] || paraD["tolerance_for_convergence"]<=0.0)
		{
			std::cerr <<"Error no or non-positive tolerance_for_convergence given for iteration_stopping_criteria=1" <<std::endl;
			return false;
		}
		break;
		case 2: break;
		default:
		std::cerr <<"Error, iteration_stopping_criteria=" <<paraI["iteration_stopping_criteria"] <<" not implemented" <<std::endl;
		return false;
	}
	if(!paraIsSet["tolerance_for_line_minimization"] || paraD["tolerance_for_line_minimization"] < 0.0)
	{
		std::cerr <<"Error, no or negative tolerance_for_line_minimization given" <<std::endl;
		return false;
	}
	if(paraIsSet["max_numer_of_iterations"] && paraI["max_numer_of_iterations"]<=0)
	{
		std::cerr <<"Error, non-positive max_numer_of_iterations given" <<std::endl;
		return false;
	}
	
	//iteration start
	if( !paraIsSet["starting_procedure"] )
	{
		std::cerr <<"Error, no starting_procedure set" <<std::endl;
		return false;
	}
	switch( paraI["starting_procedure"] )
	{
		case 1: //the same as 2
		case 2:
		if( !( paraIsSet["start_m"] && paraIsSet["start_s"] && paraIsSet["initial_step_size"] ) )
		{
			std::cerr <<"Error, no start_m, start_s or initial_step_size set" <<std::endl;
			return false;
		}
		if(paraD["initial_step_size"]<=0.0)
		{
			std::cerr <<"Error, non-positive initial_step_size chosen" <<std::endl;
			return false;
		}
		break;
		case 3://same as 4
		case 4:
		if( !( paraIsSet["minimal_m"] && paraIsSet["maximal_m"] && paraIsSet["step_m"] ) ||
		    !( paraIsSet["minimal_s"] && paraIsSet["maximal_s"] && paraIsSet["step_s"] ) )
		{
			std::cerr <<"Error, no scanning range set for initial scan ([minimal|maximal|step]_[m|s])" <<std::endl;
			return false;
		}
		if( (paraD["minimal_m"] > paraD["maximal_m"] || paraD["step_m"]<=0.0) ||
		    (paraD["minimal_s"] > paraD["maximal_s"] || paraD["step_s"]<=0.0) )
		{
			std::cerr <<"Error, inconstistent scanning range for the initial scan" <<std::endl;
			return false;
		}
		break;
		default:
		std::cerr <<"Error, starting_procedure=" <<paraI["starting_procedure"] <<" not implemented" <<std::endl;
		return false;
	}
	return true;
}//checkConsistencyOfParameters




void CEPscan_helper::fillSetWithRange( const double min, const double max, const double step, std::set< double > &toFill)
{
	toFill.clear();
	
	int numberOfEntries=static_cast< int >( (max-min)/step + 1.5);
	
	for( int i=0; i<numberOfEntries; ++i)
	{
		toFill.insert( min + static_cast< double >(i) *step);
	}
}




void CEPscan_helper::findStartByTwo_1d_Scans(std::set< double > &testValues_magnetization, std::set< double > &testValues_staggeredMagnetization, constrainedEffectivePotential &CEP, double &start_m, double &start_s)
{
	if(testValues_magnetization.empty() || testValues_staggeredMagnetization.empty())
	{
		std::cerr <<"Error, empty sets in CEPscan_helper::findStartByTwo_1d_Scans()" <<std::endl;
		exit(EXIT_FAILURE);
	}
	double smallesResult=CEP.computeConstrainedEffectivePotential_onlyFunction(*testValues_magnetization.begin(), *testValues_staggeredMagnetization.begin());
// 	std::cout <<"scan: U(" <<*testValues_magnetization.begin()<<", " <<*testValues_staggeredMagnetization.begin() <<")=" << smallesResult <<std::endl;
	start_m=*testValues_magnetization.begin();
	start_s=*testValues_staggeredMagnetization.begin();
	double testResult(0.0);
	
	for(std::set< double >::iterator magIter=++testValues_magnetization.begin(); magIter!=testValues_magnetization.end(); ++magIter)
	{
		testResult=CEP.computeConstrainedEffectivePotential_onlyFunction(*magIter, *testValues_staggeredMagnetization.begin());
// 		std::cout <<"scan: U(" <<*magIter<<", " <<*testValues_staggeredMagnetization.begin() <<")=" <<testResult <<std::endl;
		if(testResult < smallesResult)
		{
			smallesResult=testResult; 
			start_m=*magIter;
			start_s=*testValues_staggeredMagnetization.begin();
		}
	}
	for(std::set< double >::iterator stagIter=++testValues_staggeredMagnetization.begin(); stagIter!=testValues_staggeredMagnetization.end(); ++stagIter)
	{
		testResult=CEP.computeConstrainedEffectivePotential_onlyFunction(*testValues_magnetization.begin(), *stagIter);
// 		std::cout <<"scan: U(" <<*testValues_magnetization.begin()<<", " <<*stagIter <<")=" <<testResult <<std::endl;
		if(testResult < smallesResult)
		{
			smallesResult=testResult; 
			start_m=*testValues_magnetization.begin();
			start_s=*stagIter;
		}
	}
}//findStartByTwo_1d_Scans




void CEPscan_helper::findStartByOne_2d_Scan(std::set< double > &testValues_magnetization, std::set< double > &testValues_staggeredMagnetization, constrainedEffectivePotential &CEP, double &start_m, double &start_s)
{
	if(testValues_magnetization.empty() || testValues_staggeredMagnetization.empty())
	{
		std::cerr <<"Error, empty sets in CEPscan_helper::findStartByTwo_1d_Scans()" <<std::endl;
		exit(EXIT_FAILURE);
	}
	double smallesResult=CEP.computeConstrainedEffectivePotential_onlyFunction(*testValues_magnetization.begin(), *testValues_staggeredMagnetization.begin());
// 	std::cout <<"scan: U(" <<*testValues_magnetization.begin()<<", " <<*testValues_staggeredMagnetization.begin() <<")=" << smallesResult <<std::endl;
	start_m=*testValues_magnetization.begin();
	start_s=*testValues_staggeredMagnetization.begin();
	double testResult(0.0);
	
	for(std::set< double >::iterator magIter=++testValues_magnetization.begin(); magIter!=testValues_magnetization.end(); ++magIter)
	{
		for(std::set< double >::iterator stagIter=++testValues_staggeredMagnetization.begin(); stagIter!=testValues_staggeredMagnetization.end(); ++stagIter)
		{
			testResult=CEP.computeConstrainedEffectivePotential_onlyFunction(*magIter, *stagIter);
// 			std::cout <<"scan: U(" <<*magIter<<", " <<*stagIter <<")=" <<testResult <<std::endl;
			if(testResult < smallesResult)
			{
				smallesResult=testResult; 
				start_m=*magIter;
				start_s=*stagIter;
			}
		}
	}
	
}


bool CEPscan_helper::printResultsVectorToStream(const std::vector< resultForOutput > &results, std::ostream &output)
{
	//increase precision, but store the old one
	std::streamsize oldPrec=output.precision();
	output.precision(12);
	for(std::vector< resultForOutput >::const_iterator iter=results.begin(); iter!=results.end(); ++iter)
	{
		if( !( output <<iter->kappa_N <<" " <<iter->lambda_N <<" " <<iter->yukawa_N <<" " <<iter->magnetization <<" " <<iter->staggered_magnetization <<" " <<iter->potential <<" " <<iter->d2U_ov_dmdm <<" " <<iter->d2U_ov_dsds <<" " <<iter->d2U_ov_dmds  <<std::endl ) )
		{
			std::cerr <<"Error during output of results" <<std::endl;
			return false;
		}
	}
	output.precision(oldPrec);
	return true;

}





