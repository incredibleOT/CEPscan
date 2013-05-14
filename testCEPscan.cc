#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "constrainedEffectivePotential.h"


using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
	
// 	int L0(4), L1(4), L2(4), L3(4);
// 	int L0(16), L1(16), L2(16), L3(32);
// 	int L0(32), L1(32), L2(32), L3(64);
	int L0(64), L1(64), L2(64), L3(128);
// 	int L0(8), L1(10), L2(14), L3(22);
// 	bool antiL3(false);
	bool antiL3(false);
	
// 	double y=175.0/246.0;
	double y=4.0;
	
	double kappa=-0.5;
	double lambda=0.3;
	
	constrainedEffectivePotential CEP(L0,L1,L2,L3,antiL3);
	cout <<"constrainedEffectivePotential initialized" <<endl;
	CEP.set_yukawa_N(y);
	CEP.set_kappa_N(kappa);
	CEP.set_lambda_N(lambda);
	
	
	
	
	/*
	double mCentral=0.2;
	double sCentral=0.4;
	double mStep=0.0001;
	double sStep=0.0001;
	int steps=5;
	int maxCounter(2*steps+1);
	
	double *resultsMvar = new double [maxCounter];
	double *resultsSvar = new double [maxCounter];
	double m,s;
	
	
	
	//test derivative
	for(int im=0; im<maxCounter; ++im)
	{
		m=mCentral-mStep*(steps-im);
		resultsMvar[im]=CEP.computeConstrainedEffectivePotential_onlyFunction(m,sCentral);
		s=sCentral-sStep*(steps-im);
		resultsSvar[im]=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,s);
	}
	
	for(int i=0; i<maxCounter; ++i)
	{
		cout.precision(15);
		m=mCentral-mStep*(steps-i);
		s=sCentral-sStep*(steps-i);
		cout <<"U(m=" <<m <<", s=" <<sCentral <<") = " <<resultsMvar[i] <<"   U(m=" <<mCentral <<", s=" <<s <<") = " <<resultsSvar[i] <<endl;
	}
	for(int i=0; i<steps; ++i)
	{
		cout.precision(15);
		double dist_m=mStep*(i+1);
		double dist_s=sStep*(i+1);
		cout <<"(U(m=" <<mCentral <<"+" <<dist_m <<", s=" <<sCentral <<") -  U(m=" <<mCentral <<"-" <<dist_m <<", s=" <<sCentral <<"))/"<<2*dist_m <<"  = " 
		     <<(resultsMvar[steps+1+i]-resultsMvar[steps-i-1])/(2.0*dist_m) <<endl;
		     
		cout <<"(U(m=" <<mCentral <<", s=" <<sCentral <<"+" <<dist_s <<") -  U(m=" <<mCentral <<", s=" <<sCentral <<"-" <<dist_s <<"))/"<<2*dist_s <<"  = " 
		     <<(resultsSvar[steps+1+i]-resultsSvar[steps-1-i])/(2.0*dist_s) <<endl;
	}
	
	double res_U_dm(0.0), res_U_ds(0.0);
	CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral,sCentral, res_U_dm, res_U_ds);
	cout <<"dU/dm(m=" <<mCentral <<", s=" <<sCentral <<") = " <<res_U_dm <<"   dU/ds(m=" <<mCentral <<", s=" <<sCentral <<") = " <<res_U_ds <<endl;
	
	delete [] resultsMvar;
	delete [] resultsSvar;
	
	//now use the function that computes teh function and the derivative at the same time
	
	double U(0.0),U_dm(0.0),U_ds(0.0);
	for(int im=0; im<maxCounter; ++im)
	{
		m=mCentral-mStep*(steps-im);
		CEP.computeConstrainedEffectivePotential_FunctionAndGradient(m,sCentral, U, U_dm, U_ds);
		cout <<"m="<<m <<"  s=" <<sCentral <<"    U(m,s) = " <<U <<"   dU/dm = " <<U_dm <<"   dU/ds = " <<U_ds <<endl;
		s=sCentral-sStep*(steps-im);
		CEP.computeConstrainedEffectivePotential_FunctionAndGradient(mCentral,s, U, U_dm, U_ds);
		cout <<"m="<<mCentral <<"  s=" <<s <<"    U(m,s) = " <<U <<"   dU/dm = " <<U_dm <<"   dU/ds = " <<U_ds <<endl;
	}
	
	
	cout <<endl <<"Testing the gsl version of only the function:" <<endl;
	double resultNormal=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral);
	gsl_vector *testVector = gsl_vector_alloc(2);
	gsl_vector_set(testVector, 0,  mCentral);
	gsl_vector_set(testVector, 1,  sCentral);
	void *dummyPTR;
	double result_gsl=CEP.computeConstrainedEffectivePotential_onlyFunction_gsl(testVector, dummyPTR);
	cout <<"Normal result: " <<resultNormal <<"   from gsl-function: " <<result_gsl <<endl;
	
	cout <<endl <<"Testing the gsl version of only the gradient:" <<endl;
	gsl_vector *testGradient = gsl_vector_alloc(2);
	double grad_dm(0.0), grad_ds(0.0);
	CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral,sCentral, grad_dm, grad_ds);
	CEP.computeConstrainedEffectivePotential_onlyGradient_gsl(testVector, dummyPTR, testGradient);
	cout <<"dU/dm normal: " <<grad_dm <<"   from gsl: " <<gsl_vector_get(testGradient,0) <<endl;
	cout <<"dU/ds normal: " <<grad_ds <<"   from gsl: " <<gsl_vector_get(testGradient,1) <<endl;
	
	cout <<endl <<"Testing the gsl version of function and gradient:" <<endl;
	
	grad_dm=0.0, grad_ds=0.0, resultNormal=0.0;
	result_gsl=0;gsl_vector_set_zero(testGradient);
	CEP.computeConstrainedEffectivePotential_FunctionAndGradient(mCentral,sCentral, resultNormal, grad_dm, grad_ds);
	CEP.computeConstrainedEffectivePotential_FunctionAndGradient_gsl(testVector, dummyPTR, &result_gsl, testGradient);
	cout <<"Normal result: " <<resultNormal <<"   from gsl-function: " <<result_gsl <<endl;
	cout <<"dU/dm normal: " <<grad_dm <<"   from gsl: " <<gsl_vector_get(testGradient,0) <<endl;
	cout <<"dU/ds normal: " <<grad_ds <<"   from gsl: " <<gsl_vector_get(testGradient,1) <<endl;
	
	
	
	cout <<endl <<"Testing the gsl wrapper of function:" <<endl;
	grad_dm=0.0, grad_ds=0.0, resultNormal=0.0;
	result_gsl=0;gsl_vector_set_zero(testGradient);
	
	constrainedEffectivePotential *CEP_ptr=&CEP;
	cout <<"wrapper returns: " <<wrapper_computeConstrainedEffectivePotential_onlyFunction_gsl(testVector, (void *)CEP_ptr) <<endl;
	
	wrapper_computeConstrainedEffectivePotential_onlyGradient_gsl(testVector, CEP_ptr, testGradient);
	cout <<"dU/dm from wrapper_gsl_onlyGradient: " <<gsl_vector_get(testGradient,0) <<endl;
	cout <<"dU/ds from wrapper_gsl_onlyGradient: " <<gsl_vector_get(testGradient,1) <<endl;
	
	grad_dm=0.0, grad_ds=0.0, resultNormal=0.0;
	result_gsl=0;gsl_vector_set_zero(testGradient);
	
	wrapper_computeConstrainedEffectivePotential_FunctionAndGradient_gsl(testVector, CEP_ptr, &result_gsl, testGradient);
	cout <<"result from wrapper_gsl_functionAndGradient: " <<result_gsl <<endl;
	cout <<"dU/dm from wrapper_gsl_functionAndGradient: " <<gsl_vector_get(testGradient,0) <<endl;
	cout <<"dU/ds from wrapper_gsl_functionAndGradient: " <<gsl_vector_get(testGradient,1) <<endl;
	*/
	
	
	
	double toleranceForLine=0.00001;
	double toleranceForConv=0.001;
	double initStepSize=0.01;
	int miniAlgo=3;
	CEP.set_toleranceForLineMinimization(toleranceForLine);
	CEP.set_toleranceForConvergence(toleranceForConv);
	CEP.set_initialStepSize(initStepSize);
	CEP.set_minimizationAlgorithm(miniAlgo);
	
	int returnFromInit=CEP.initializeMinimizer(1.1,10.);
	cout <<"initialization returned: " <<returnFromInit <<endl;
	
// 	test iterator
	double actMag(0.0),actStagMag(0.0),actPot(0.0),act_dU_ov_dm(0.0), act_dU_ov_ds(0.0);
	CEP.reInitializeMinimizer();

	
	CEP.getActualMinimizerLocation(actMag, actStagMag),
	actPot=CEP.getActualMinimizerValue();
	CEP.getActualMinimizerGradient(act_dU_ov_dm, act_dU_ov_ds);
	
	cout <<"After initialization:" <<endl;
	cout <<"mag=" <<actMag <<"  stag=" <<actStagMag <<"  U=" <<actPot <<"  dU/dm=" <<act_dU_ov_dm <<"  dU/ds=" <<act_dU_ov_ds <<endl;
	
	for(int i=0; i<100; ++i)
	{
		CEP.iterateMinimizer();
		CEP.getActualMinimizerLocation(actMag, actStagMag);
		actPot=CEP.getActualMinimizerValue();
		CEP.getActualMinimizerGradient(act_dU_ov_dm, act_dU_ov_ds);
		cout <<"Iteration number " <<i <<", iterationStopped() returns: " <<CEP.iterationStopped() <<", testMinimizerGradient() returns: " <<CEP.testMinimizerGradient() <<endl;
		cout <<"mag=" <<actMag <<"  stag=" <<actStagMag <<"  U=" <<actPot <<"  dU/dm=" <<act_dU_ov_dm <<"  dU/ds=" <<act_dU_ov_ds <<endl;
		if(CEP.iterationStopped()){break;}
	}
	
	CEP.reInitializeMinimizer();

	CEP.getActualMinimizerLocation(actMag, actStagMag),
	actPot=CEP.getActualMinimizerValue();
	CEP.getActualMinimizerGradient(act_dU_ov_dm, act_dU_ov_ds);
	
	cout <<"After reInitialization:" <<endl;
	cout <<"mag=" <<actMag <<"  stag=" <<actStagMag <<"  U=" <<actPot <<"  dU/dm=" <<act_dU_ov_dm <<"  dU/ds=" <<act_dU_ov_ds <<endl;
	
	
	double k_min=-0.5, k_max=0.2, k_step=0.001;
	double actualKappa=k_min;
	/*
	std::map< double, std::pair< double, double > > minimum_fromLast; //takes last result
	std::map< double, std::pair< double, double > > minimum_fromMag; //starts at 2.0,0.001
	std::map< double, std::pair< double, double > > minimum_fromStag; //starts at 0.001,2.0
	std::map< double, double > function_fromLast;
	std::map< double, double > function_fromMag;
	std::map< double, double > function_fromStag;
// 	std::map< double, std::pair< double, double > > decreaseKappaRun;
	//forward
	cout <<"measurement for increasing kappa" <<endl;
	cout <<"start at last value" <<endl;
	CEP.reInitializeMinimizer(1.0,1.0);
	actualKappa=k_min;
	while(actualKappa<=k_max)
	{
		CEP.set_kappa_N(actualKappa);
		int retFromIter=CEP.iterateUntilIterationStopps();
		if(retFromIter)
		{
			double m,s;
			CEP.getActualMinimizerLocation(m,s);
// 			cout <<"kappa=" <<CEP.get_kappa_N() <<"   m=" <<m <<"   s=" <<s <<",   after " <<retFromIter <<" iterations" <<endl;
			minimum_fromLast.insert( std::make_pair(CEP.get_kappa_N(), std::make_pair(m,s)));
			function_fromLast.insert(std::make_pair(CEP.get_kappa_N(), CEP.getActualMinimizerValue() ) );
		}
		actualKappa+=k_step;
	}
	cout <<"start at 2.0, 0.001" <<endl;
	actualKappa=k_min;
	while(actualKappa<=k_max)
	{
		CEP.set_kappa_N(actualKappa);
		CEP.reInitializeMinimizer(2.0,0.001);
		int retFromIter=CEP.iterateUntilIterationStopps();
		if(retFromIter)
		{
			double m,s;
			CEP.getActualMinimizerLocation(m,s);
// 			cout <<"kappa=" <<CEP.get_kappa_N() <<"   m=" <<m <<"   s=" <<s <<",   after " <<retFromIter <<" iterations" <<endl;
			minimum_fromMag.insert( std::make_pair(CEP.get_kappa_N(), std::make_pair(m,s)));
			function_fromMag.insert(std::make_pair(CEP.get_kappa_N(), CEP.getActualMinimizerValue() ) );
		}
		else{ cout <<"No convergence for kappa=" <<CEP.get_kappa_N() <<endl; }
		actualKappa+=k_step;
	}
	cout <<"start at 0.001, 2.0" <<endl;
	actualKappa=k_min;
	while(actualKappa<=k_max)
	{
		CEP.set_kappa_N(actualKappa);
		CEP.reInitializeMinimizer(0.001,2.0);
		int retFromIter=CEP.iterateUntilIterationStopps();
		if(retFromIter)
		{
			double m,s;
			CEP.getActualMinimizerLocation(m,s);
// 			cout <<"kappa=" <<CEP.get_kappa_N() <<"   m=" <<m <<"   s=" <<s <<",   after " <<retFromIter <<" iterations" <<endl;
			minimum_fromStag.insert( std::make_pair(CEP.get_kappa_N(), std::make_pair(m,s)));
			function_fromStag.insert(std::make_pair(CEP.get_kappa_N(), CEP.getActualMinimizerValue() ) );
		}
		actualKappa+=k_step;
	}
	
	//fromLast
	{
		std::ofstream outputFile( "output/testoutput_fromLast_3.txt" );
		if(!outputFile.good())
		{
			cerr <<"Error opening output file" <<endl;
			exit(EXIT_FAILURE);
		}
		outputFile.precision(15);
		for( std::map< double, std::pair< double, double > >::const_iterator iter=minimum_fromLast.begin(); iter!=minimum_fromLast.end(); ++iter)
		{
			outputFile <<iter->first <<"  " <<iter->second.first <<"  " <<iter->second.second <<"  " <<function_fromLast[iter->first] <<endl;
		}
		outputFile.close();
	}
	//fromMag
	{
		std::ofstream outputFile( "output/testoutput_fromMag_3.txt" );
		if(!outputFile.good())
		{
			cerr <<"Error opening output file" <<endl;
			exit(EXIT_FAILURE);
		}
		outputFile.precision(15);
		for( std::map< double, std::pair< double, double > >::const_iterator iter=minimum_fromMag.begin(); iter!=minimum_fromMag.end(); ++iter)
		{
			outputFile <<iter->first <<"  " <<iter->second.first <<"  " <<iter->second.second <<"  " <<function_fromMag[iter->first] <<endl;
		}
		outputFile.close();
	}
	//fromStag
	{
		std::ofstream outputFile( "output/testoutput_fromStag_3.txt" );
		if(!outputFile.good())
		{
			cerr <<"Error opening output file" <<endl;
			exit(EXIT_FAILURE);
		}
		outputFile.precision(15);
		for( std::map< double, std::pair< double, double > >::const_iterator iter=minimum_fromStag.begin(); iter!=minimum_fromStag.end(); ++iter)
		{
			outputFile <<iter->first <<"  " <<iter->second.first <<"  " <<iter->second.second <<"  " <<function_fromStag[iter->first] <<endl;
		}
		outputFile.close();
	}
	*/
	actualKappa=-0.383;
	CEP.set_kappa_N(actualKappa);
	{
		std::ofstream outputFile( "output/testoutput_2d_scan_k-0.383_L64.txt" );
		if(!outputFile.good())
		{
			cerr <<"Error opening output file" <<endl;
			exit(EXIT_FAILURE);
		}
		outputFile.precision(15);
		double m_min=0.006, m_max=0.006, m_step=0.0002;
		double s_min=1.7, s_max=2.05, s_step=0.0002;
		
		double actual_m=m_min;
		while(actual_m<=m_max)
		{
			cout <<"evaluating mag=" <<actual_m <<endl;
			double actual_s=s_min;
			while(actual_s<=s_max)
			{
				double result=CEP.computeConstrainedEffectivePotential_onlyFunction(actual_m, actual_s);
				outputFile <<actual_m <<"  " <<actual_s <<"  " <<result <<endl;
				actual_s+=s_step;
			}
			actual_m+=m_step;
		}
		outputFile.close();
	}

	

}
