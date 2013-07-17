#include "plotPotential_helper.h"



void plotPotential_helper::prepareParameterMaps( std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, std::string > &paraS, std::map< std::string, bool > &paraIsSet )
{
	paraD.clear(); paraI.clear(); paraS.clear(); paraIsSet.clear();
	
	paraI["L_s"]=-1;
	paraI["L_t"]=-1;
	paraI["antiperiodic_L_t"]=-1;
	
	paraI["include_bosonic_loop"]=0;
	paraI["use_improved_treeLevel"]=0;
	paraI["use_improved_firstOrder"]=0;
	
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
	
	paraI["scan_lambda_6_N"]=0;
	paraD["lambda_6_N"] = 0.0;
	paraD["lambda_6_N_min"] = 0.0;
	paraD["lambda_6_N_max"] = 0.0;
	paraD["lambda_6_N_step"] = 0.0;
	
	paraI["scan_yukawa_N"]=0;
	paraD["yukawa_N"] = 0.0;
	paraD["yukawa_N_min"] = 0.0;
	paraD["yukawa_N_max"] = 0.0;
	paraD["yukawa_N_step"] = 0.0;
	
	paraI["scan_mag"]=0;
	paraD["mag"] = 0.0;
	paraD["mag_min"] = 0.0;
	paraD["mag_max"] = 0.0;
	paraD["mag_step"] = 0.0;
	
	paraI["scan_stag"]=0;
	paraD["stag"] = 0.0;
	paraD["stag_min"] = 0.0;
	paraD["stag_max"] = 0.0;
	paraD["stag_step"] = 0.0;
	
	paraI["interpret_yukawa_as_continuum"]=0;
	paraI["interpret_lambda_as_continuum"]=0;
	paraI["interpret_lambda_6_as_continuum"]=0;
	
	//directly set to default values
	paraI["N_f"]=1;
	paraD["rho"]=1.0;
	paraD["r"]=0.5;
	
	
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


std::string plotPotential_helper::generateOutputFileName( const std::string &original, std::map< std::string, double > &paraD, std::map< std::string, int > &paraI, std::map< std::string, bool > &paraIsSet )
{
	std::string outputFileName(original);
	if( outputFileName.find("[Ls]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"L" <<paraI["L_s"];
		outputFileName.replace(outputFileName.find("[Ls]"),4, ss.str() );
	}
	if( outputFileName.find("[Lt]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"T" <<paraI["L_t"];
		outputFileName.replace(outputFileName.find("[Lt]"),4, ss.str() );
	}
	if( outputFileName.find("[k]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"k_";
		ss.precision(7);
		if(paraI["scan_kappa_N"]){ss <<paraD["kappa_N_min"] <<"_"<<paraD["kappa_N_max"]; }
		else {ss <<paraD["kappa_N"]; }
		outputFileName.replace(outputFileName.find("[k]"),3, ss.str() );
	}
	if( outputFileName.find("[l]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"l_";
		ss.precision(7);
		if(paraI["scan_lambda_N"]){ss <<paraD["lambda_N_min"] <<"_"<<paraD["lambda_N_max"]; }
		else {ss <<paraD["lambda_N"]; }
		outputFileName.replace(outputFileName.find("[l]"),3, ss.str() );
	}
	if( outputFileName.find("[l6]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"l6_";
		ss.precision(7);
		if(paraI["scan_lambda_6_N"]){ss <<paraD["lambda_6_N_min"] <<"_"<<paraD["lambda_6_N_max"]; }
		else {ss <<paraD["lambda_6_N"]; }
		outputFileName.replace(outputFileName.find("[l6]"),4, ss.str() );
	}
	if( outputFileName.find("[y]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"y_";
		ss.precision(7);
		if(paraI["scan_yukawa_N"]){ss <<paraD["yukawa_N_min"] <<"_"<<paraD["yukawa_N_max"]; }
		else {ss <<paraD["yukawa_N"]; }
		outputFileName.replace(outputFileName.find("[y]"),3, ss.str() );
	}
	if( outputFileName.find("[m]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"mag_";
		ss.precision(7);
		if(paraI["scan_mag"]){ss <<paraD["mag_min"] <<"_"<<paraD["mag_max"]; }
		else {ss <<paraD["mag"]; }
		outputFileName.replace(outputFileName.find("[m]"),3, ss.str() );
	}
	if( outputFileName.find("[s]")!=std::string::npos )
	{
		std::ostringstream ss;
		ss <<"stag_";
		ss.precision(7);
		if(paraI["scan_stag"]){ss <<paraD["stag_min"] <<"_"<<paraD["stag_max"]; }
		else {ss <<paraD["stag"]; }
		outputFileName.replace(outputFileName.find("[s]"),3, ss.str() );
	}
	if( outputFileName.find("[loop]")!=std::string::npos )
	{
		std::ostringstream ss;
		if(paraI["include_bosonic_loop"]){ ss <<"withLoop";}
		else{ ss <<"noLoop";}
		outputFileName.replace(outputFileName.find("[loop]"),6, ss.str() );
	}
	if( outputFileName.find("[bc]")!=std::string::npos )
	{
		std::ostringstream ss;
		if(paraI["antiperiodic_L_t"]){ ss <<"aBC";}
		else{ ss <<"pBC";}
		outputFileName.replace(outputFileName.find("[bc]"),4, ss.str() );
	}
	if( outputFileName.find("[det]")!=std::string::npos )
	{
		std::ostringstream ss;
		if(paraI["use_improved_treeLevel"]){ ss <<"withDet";}
		else{ ss <<"noDet";}
		outputFileName.replace(outputFileName.find("[det]"),5, ss.str() );
	}
	if( outputFileName.find("[1st]")!=std::string::npos )
	{
		std::ostringstream ss;
		if(paraI["use_improved_firstOrder"]){ ss <<"with1st";}
		else{ ss <<"no1st";}
		outputFileName.replace(outputFileName.find("[1st]"),5, ss.str() );
	}
	
	return outputFileName;
}//generateOutputFileName


bool plotPotential_helper::printResultsVectorToStream(const std::vector< resultForOutput > &results, std::ostream &output)
{
	//increase precision, but store the old one
	std::streamsize oldPrec=output.precision();
	output.precision(12);
	for(std::vector< resultForOutput >::const_iterator iter=results.begin(); iter!=results.end(); ++iter)
	{
		if( !( output <<iter->kappa_N <<" " <<iter->lambda_N <<" " <<iter->lambda_6_N <<" " <<iter->yukawa_N <<" " <<iter->magnetization <<" " <<iter->staggered_magnetization <<" " <<iter->potential <<" " <<iter->dU_ov_dm <<" " <<iter->dU_ov_ds <<" " <<iter->d2U_ov_dmdm <<" " <<iter->d2U_ov_dsds <<" " <<iter->d2U_ov_dmds  <<std::endl ) )
		{
			std::cerr <<"Error during output of results" <<std::endl;
			return false;
		}
	}
	output.precision(oldPrec);
	return true;

}
