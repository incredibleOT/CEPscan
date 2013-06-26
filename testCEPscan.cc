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
	int L0(16), L1(18), L2(20), L3(32);
// 	int L0(32), L1(32), L2(32), L3(64);
// 	int L0(64), L1(64), L2(64), L3(128);
// 	int L0(8), L1(10), L2(14), L3(22);
// 	bool antiL3(false);
	bool antiL3(false);
	
// 	double y=175.0/246.0;
	double y=1.0;
	
	double kappa=+0.1;
	double lambda=0.2;
	double lambda_6=0.2;
	
	constrainedEffectivePotential CEP(L0,L1,L2,L3,antiL3);
	cout <<"constrainedEffectivePotential initialized" <<endl;
	CEP.set_yukawa_N(y);
	CEP.set_kappa_N(kappa);
	CEP.set_lambda_N(lambda);
	CEP.set_lambda_6_N(lambda_6);
	CEP.set_useBosonicLoop(false);
	
	CEP.set_useImprovedGaussian(true);
	
	//test bosonic determinant
	if(true)
	{
		double mCentral=0.2;
		double sCentral=0.4;
		int NumberOfIterations=1;
		
		double qad_determinant(0);
		for(int i=0; i<NumberOfIterations; ++i)
		{
			qad_determinant=CEP.computeBosonicDeterminantContributionForImprovedGaussian_onlyFunction_qad(mCentral, sCentral);
		}
		cout <<endl <<"Bosonic derterminant qad: " <<qad_determinant <<endl;
		
		double fromStored_determinant(0);
		for(int i=0; i<NumberOfIterations; ++i)
		{
			fromStored_determinant=CEP.computeBosonicDeterminantContributionForImprovedGaussian_onlyFunction_fromStoredSumOfCos(mCentral, sCentral);
		}
		cout <<endl <<"Bosonic derterminant from stored: " <<fromStored_determinant <<endl;
		
		double fromStored_dDet_dm(0.0), fromStored_dDet_ds(0.0);
		for(int i=0; i<NumberOfIterations; ++i)
		{
			CEP.computeBosonicDeterminantContributionForImprovedGaussian_onlyGradient_fromStoredSumOfCos(mCentral, sCentral, fromStored_dDet_dm, fromStored_dDet_ds);
		}
		cout <<endl <<"Bosonic derterminant_ov_dm from stored: " <<fromStored_dDet_dm <<"    Bosonic derterminant_ov_ds from stored: " <<fromStored_dDet_ds <<endl;
		
		double fromStoredComb_U, fromStoredComb_dm, fromStoredComb_ds;
		CEP.computeBosonicDeterminantContributionForImprovedGaussian_FunctionAndGradient_fromStoredSumOfCos(mCentral, sCentral, fromStoredComb_U, fromStoredComb_dm, fromStoredComb_ds);
		cout <<endl <<"bosDet from comb: " <<fromStoredComb_U <<"   dm:" <<fromStoredComb_dm <<"    ds: " <<fromStoredComb_ds <<endl;
		
		
	}
	
	
	
	
	//test derivatives
	if(true)
	{
		cout.precision(15);
		double mCentral(0.3), sCentral(0.4);
		double eps=0.00001;
		double result(0.0), dUdm(0.0), dUds(0.0), dummy(0.0);
		double dUdmdm(0.0), dUdsds(0.0), dUdmds(0.0);
		double eps1(0.00001), eps2(0.000001), eps3(0.0000001);
		double dUdmdm1(0.0), dUdsds1(0.0), dUdmds1(0.0);
		double dUdmdm2(0.0), dUdsds2(0.0), dUdmds2(0.0);
		double dUdmdm3(0.0), dUdsds3(0.0), dUdmds3(0.0);
		
		double dUdm1(0.0), dUds1(0.0);
		double dUdm2(0.0), dUds2(0.0);
		double dUdm3(0.0), dUds3(0.0);
	
		cout <<"From function only:" <<endl;
		result=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral);
		cout <<"U = " << result <<endl;
		dUdm=(CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral+eps,sCentral) - CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral-eps,sCentral))/(2.0*eps);
		cout <<"dU/dm= " <<dUdm <<endl;
		dUds=(CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral+eps) - CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral-eps))/(2.0*eps);
		cout <<"dU/ds= " <<dUds <<endl;
		
		cout <<endl <<"From gradient only:" <<endl;
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral, sCentral, dUdm, dUds);
		cout <<"dU/dm= " <<dUdm <<endl;
		cout <<"dU/ds= " <<dUds <<endl;
		
		cout <<endl <<"From FunctionAndGradient:" <<endl;
		CEP.computeConstrainedEffectivePotential_FunctionAndGradient(mCentral, sCentral, result, dUdm, dUds);
		cout <<"U = " << result <<endl;
		cout <<"dU/dm= " <<dUdm <<endl;
		cout <<"dU/ds= " <<dUds <<endl;
		
		
		cout <<endl <<"First derivatives from the function" <<endl;
		
		
		result=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral+eps1,sCentral);
		dummy=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral-eps1,sCentral);
		dUdm1=0.5/eps1*(result-dummy);
		result=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral+eps1);
		dummy=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral-eps1);
		dUds1=0.5/eps1*(result-dummy);
		
		result=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral+eps2,sCentral);
		dummy=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral-eps2,sCentral);
		dUdm2=0.5/eps2*(result-dummy);
		result=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral+eps2);
		dummy=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral-eps2);
		dUds2=0.5/eps2*(result-dummy);
		
		result=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral+eps3,sCentral);
		dummy=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral-eps3,sCentral);
		dUdm3=0.5/eps3*(result-dummy);
		result=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral+eps3);
		dummy=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral-eps3);
		dUds3=0.5/eps3*(result-dummy);
		
		cout <<"eps=" <<eps1 <<"   dU/dm=" <<dUdm1 <<"   dU/ds=" <<dUds1 <<endl;
		cout <<"eps=" <<eps2 <<"   dU/dm=" <<dUdm2 <<"   dU/ds=" <<dUds2 <<endl;
		cout <<"eps=" <<eps3 <<"   dU/dm=" <<dUdm3 <<"   dU/ds=" <<dUds3 <<endl;
		
		
		cout <<endl <<"First derivatives from the gradient" <<endl;
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral,sCentral,dUdm,dUds);
		cout <<"dU/dm=" <<dUdm <<"   dU/ds=" <<dUds <<endl;
		
		cout <<endl <<"relative difference" <<endl;
		
		cout <<"eps=" <<eps1 <<"   dU/dm:" <<0.5*(dUdm1-dUdm)/(dUdm1+dUdm);
		cout <<"   dU/ds:" <<0.5*(dUds1-dUds)/(dUds1+dUds) <<endl;
		
		cout <<"eps=" <<eps2 <<"   dU/dm:" <<0.5*(dUdm2-dUdm)/(dUdm2+dUdm);
		cout <<"   dU/ds:" <<0.5*(dUds2-dUds)/(dUds2+dUds) <<endl;
		
		cout <<"eps=" <<eps3 <<"   dU/dm:" <<0.5*(dUdm3-dUdm)/(dUdm3+dUdm);
		cout <<"   dU/ds:" <<0.5*(dUds3-dUds)/(dUds3+dUds) <<endl;
		
	
		
		
		
		cout <<endl <<"Second derivatives from the gradient" <<endl;
		
		
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral+eps1, sCentral, dUdm, dUds);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral-eps1, sCentral, dUdm2, dUds2);
		dUdmdm1=0.5/eps1*(dUdm-dUdm2);
		dUdmds1=0.5/eps1*(dUds-dUds2);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral, sCentral+eps1, dUdm, dUds);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral, sCentral-eps1, dUdm2, dUds2);
		dUdsds1=0.5/eps1*(dUds-dUds2);
		
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral+eps2, sCentral, dUdm, dUds);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral-eps2, sCentral, dUdm2, dUds2);
		dUdmdm2=0.5/eps2*(dUdm-dUdm2);
		dUdmds2=0.5/eps2*(dUds-dUds2);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral, sCentral+eps2, dUdm, dUds);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral, sCentral-eps2, dUdm2, dUds2);
		dUdsds2=0.5/eps2*(dUds-dUds2);
		
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral+eps3, sCentral, dUdm, dUds);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral-eps3, sCentral, dUdm2, dUds2);
		dUdmdm3=0.5/eps3*(dUdm-dUdm2);
		dUdmds3=0.5/eps3*(dUds-dUds2);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral, sCentral+eps3, dUdm, dUds);
		CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral, sCentral-eps3, dUdm2, dUds2);
		dUdsds3=0.5/eps3*(dUds-dUds2);
		
		cout <<"eps=" <<eps1 <<"   d/dm(dU/dm)=" <<dUdmdm1 <<"   d/ds(dU/ds)=" <<dUdsds1 <<"   d/dm(dU/ds)=" <<dUdmds1 <<endl;
		cout <<"eps=" <<eps2 <<"   d/dm(dU/dm)=" <<dUdmdm2 <<"   d/ds(dU/ds)=" <<dUdsds2 <<"   d/dm(dU/ds)=" <<dUdmds2 <<endl;
		cout <<"eps=" <<eps3 <<"   d/dm(dU/dm)=" <<dUdmdm3 <<"   d/ds(dU/ds)=" <<dUdsds3 <<"   d/dm(dU/ds)=" <<dUdmds3 <<endl;
		
		cout <<endl <<"Second derivatives from second derivative" <<endl;
		CEP.computeConstrainedEffectivePotential_secondDerivatives(mCentral, sCentral, dUdmdm, dUdsds, dUdmds);
		cout <<"d^2U/dm^2=" <<dUdmdm <<"   d^2U/ds^2=" <<dUdsds <<"   d^2U/dmds= " <<dUdmds <<endl;
		
		cout <<endl <<"relative difference" <<endl;
		
		cout <<"eps=" <<eps1 <<"   d^2U/dm^2:" <<0.5*(dUdmdm1-dUdmdm)/(dUdmdm1+dUdmdm);
		cout <<"   d^2U/dm^2:" <<0.5*(dUdsds1-dUdsds)/(dUdsds1+dUdsds);
		cout <<"   d^2U/dmds:" <<0.5*(dUdmds1-dUdmds)/(dUdmds1+dUdmds) <<endl;
		cout <<"eps=" <<eps2 <<"   d^2U/dm^2:" <<0.5*(dUdmdm2-dUdmdm)/(dUdmdm2+dUdmdm);
		cout <<"   d^2U/dm^2:" <<0.5*(dUdsds2-dUdsds)/(dUdsds2+dUdsds);
		cout <<"   d^2U/dmds:" <<0.5*(dUdmds2-dUdmds)/(dUdmds2+dUdmds) <<endl;
		cout <<"eps=" <<eps3 <<"   d^2U/dm^2:" <<0.5*(dUdmdm3-dUdmdm)/(dUdmdm3+dUdmdm);
		cout <<"   d^2U/dm^2:" <<0.5*(dUdsds3-dUdsds)/(dUdsds3+dUdsds);
		cout <<"   d^2U/dmds:" <<0.5*(dUdmds3-dUdmds)/(dUdmds3+dUdmds) <<endl;
		
		
		
	}
	
	
// 	//test derivative
// 	for(int im=0; im<maxCounter; ++im)
// 	{
// 		m=mCentral-mStep*(steps-im);
// 		resultsMvar[im]=CEP.computeConstrainedEffectivePotential_onlyFunction(m,sCentral);
// 		s=sCentral-sStep*(steps-im);
// 		resultsSvar[im]=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,s);
// 	}
// 	
// 	for(int i=0; i<maxCounter; ++i)
// 	{
// 		cout.precision(15);
// 		m=mCentral-mStep*(steps-i);
// 		s=sCentral-sStep*(steps-i);
// 		cout <<"U(m=" <<m <<", s=" <<sCentral <<") = " <<resultsMvar[i] <<"   U(m=" <<mCentral <<", s=" <<s <<") = " <<resultsSvar[i] <<endl;
// 	}
// 	for(int i=0; i<steps; ++i)
// 	{
// 		cout.precision(15);
// 		double dist_m=mStep*(i+1);
// 		double dist_s=sStep*(i+1);
// 		cout <<"(U(m=" <<mCentral <<"+" <<dist_m <<", s=" <<sCentral <<") -  U(m=" <<mCentral <<"-" <<dist_m <<", s=" <<sCentral <<"))/"<<2*dist_m <<"  = " 
// 		     <<(resultsMvar[steps+1+i]-resultsMvar[steps-i-1])/(2.0*dist_m) <<endl;
// 		     
// 		cout <<"(U(m=" <<mCentral <<", s=" <<sCentral <<"+" <<dist_s <<") -  U(m=" <<mCentral <<", s=" <<sCentral <<"-" <<dist_s <<"))/"<<2*dist_s <<"  = " 
// 		     <<(resultsSvar[steps+1+i]-resultsSvar[steps-1-i])/(2.0*dist_s) <<endl;
// 	}
	
// 	double res_U_dm(0.0), res_U_ds(0.0);
// 	CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral,sCentral, res_U_dm, res_U_ds);
// 	cout <<"dU/dm(m=" <<mCentral <<", s=" <<sCentral <<") = " <<res_U_dm <<"   dU/ds(m=" <<mCentral <<", s=" <<sCentral <<") = " <<res_U_ds <<endl;
// 	
// 	delete [] resultsMvar;
// 	delete [] resultsSvar;
// 	
// 	//now use the function that computes teh function and the derivative at the same time
// 	
// 	double U(0.0),U_dm(0.0),U_ds(0.0);
// 	for(int im=0; im<maxCounter; ++im)
// 	{
// 		m=mCentral-mStep*(steps-im);
// 		CEP.computeConstrainedEffectivePotential_FunctionAndGradient(m,sCentral, U, U_dm, U_ds);
// 		cout <<"m="<<m <<"  s=" <<sCentral <<"    U(m,s) = " <<U <<"   dU/dm = " <<U_dm <<"   dU/ds = " <<U_ds <<endl;
// 		s=sCentral-sStep*(steps-im);
// 		CEP.computeConstrainedEffectivePotential_FunctionAndGradient(mCentral,s, U, U_dm, U_ds);
// 		cout <<"m="<<mCentral <<"  s=" <<s <<"    U(m,s) = " <<U <<"   dU/dm = " <<U_dm <<"   dU/ds = " <<U_ds <<endl;
// 	}
// 	
	
// 	cout <<endl <<"Testing the gsl version of only the function:" <<endl;
// 	double resultNormal=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,sCentral);
// 	gsl_vector *testVector = gsl_vector_alloc(2);
// 	gsl_vector_set(testVector, 0,  mCentral);
// 	gsl_vector_set(testVector, 1,  sCentral);
// 	void *dummyPTR;
// 	double result_gsl=CEP.computeConstrainedEffectivePotential_onlyFunction_gsl(testVector, dummyPTR);
// 	cout <<"Normal result: " <<resultNormal <<"   from gsl-function: " <<result_gsl <<endl;
// 	
// 	cout <<endl <<"Testing the gsl version of only the gradient:" <<endl;
// 	gsl_vector *testGradient = gsl_vector_alloc(2);
// 	double grad_dm(0.0), grad_ds(0.0);
// 	CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral,sCentral, grad_dm, grad_ds);
// 	CEP.computeConstrainedEffectivePotential_onlyGradient_gsl(testVector, dummyPTR, testGradient);
// 	cout <<"dU/dm normal: " <<grad_dm <<"   from gsl: " <<gsl_vector_get(testGradient,0) <<endl;
// 	cout <<"dU/ds normal: " <<grad_ds <<"   from gsl: " <<gsl_vector_get(testGradient,1) <<endl;
// 	
// 	cout <<endl <<"Testing the gsl version of function and gradient:" <<endl;
// 	
// 	grad_dm=0.0, grad_ds=0.0, resultNormal=0.0;
// 	result_gsl=0;gsl_vector_set_zero(testGradient);
// 	CEP.computeConstrainedEffectivePotential_FunctionAndGradient(mCentral,sCentral, resultNormal, grad_dm, grad_ds);
// 	CEP.computeConstrainedEffectivePotential_FunctionAndGradient_gsl(testVector, dummyPTR, &result_gsl, testGradient);
// 	cout <<"Normal result: " <<resultNormal <<"   from gsl-function: " <<result_gsl <<endl;
// 	cout <<"dU/dm normal: " <<grad_dm <<"   from gsl: " <<gsl_vector_get(testGradient,0) <<endl;
// 	cout <<"dU/ds normal: " <<grad_ds <<"   from gsl: " <<gsl_vector_get(testGradient,1) <<endl;
// 	
// 	
// 	
// 	cout <<endl <<"Testing the gsl wrapper of function:" <<endl;
// 	grad_dm=0.0, grad_ds=0.0, resultNormal=0.0;
// 	result_gsl=0;gsl_vector_set_zero(testGradient);
// 	
// 	constrainedEffectivePotential *CEP_ptr=&CEP;
// 	cout <<"wrapper returns: " <<wrapper_computeConstrainedEffectivePotential_onlyFunction_gsl(testVector, (void *)CEP_ptr) <<endl;
// 	
// 	wrapper_computeConstrainedEffectivePotential_onlyGradient_gsl(testVector, CEP_ptr, testGradient);
// 	cout <<"dU/dm from wrapper_gsl_onlyGradient: " <<gsl_vector_get(testGradient,0) <<endl;
// 	cout <<"dU/ds from wrapper_gsl_onlyGradient: " <<gsl_vector_get(testGradient,1) <<endl;
// 	
// 	grad_dm=0.0, grad_ds=0.0, resultNormal=0.0;
// 	result_gsl=0;gsl_vector_set_zero(testGradient);
// 	
// 	wrapper_computeConstrainedEffectivePotential_FunctionAndGradient_gsl(testVector, CEP_ptr, &result_gsl, testGradient);
// 	cout <<"result from wrapper_gsl_functionAndGradient: " <<result_gsl <<endl;
// 	cout <<"dU/dm from wrapper_gsl_functionAndGradient: " <<gsl_vector_get(testGradient,0) <<endl;
// 	cout <<"dU/ds from wrapper_gsl_functionAndGradient: " <<gsl_vector_get(testGradient,1) <<endl;
// 	
	
	/*
	
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
	*/
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
	
	/*
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
	*/
	

	

}
