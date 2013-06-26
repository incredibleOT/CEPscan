#include "constrainedEffectivePotential.h"


void constrainedEffectivePotential::computeConstrainedEffectivePotential_onlyGradient( const double magnetization, const double staggeredMagnetization, double &dU_ov_dm, double &dU_ov_ds)
{
	//U= -8*kappa*(m^2-s^2) + m^2 + s^2 + lambda*( m^4 + s^4 + 6*m^2s^2 - 2*(m^2+s^2)) + lambda_6(m^6 + s^6 + 15*m^4*s^2 + 15*m^2*s^4) U_f(m,s)
	//dU/dm = dU_f/d_m - 16 kappa * m + 2*m + lambda*(4*m^3 + 12 *m*s^2 - 4*m) + lambda_6*(6*m^5 + 60*m^3*s^2 + 30*m*s^4)
	//dU/ds = dU_f/d_m + 16 kappa * s + 2*s + lambda*(4*s^3 + 12 *m^2*s - 4*s) + lambda_6*(6*s^5 + 30*m^4*s + 60*m^2*s^3)
	//
	//NOTE there are two possible improvements for the potential
	//
	//the 1st oneis a term coming from the boson loop, or better: the first order in lambda:
	//U+=16*(m^2+s^2)*lambda_N*N_f^{-1}*P_B
	//P_B= 1/V * sum_{p} 1/(2-4*lambda_N-4*kappa*sum{cos(P_mu)}) excluding zero and staggered mode
	//With that:
	//dU/dm +=32*m*lambda_N*N_f^{-1}*P_B
	//dU/ds +=32*s*lambda_N*N_f^{-1}*P_B 
	//NOTE, does not work with lambda_6!=0
	//
	//the 2nd possibility is to include an improved bosonic determinant:
	//U+=-1/(2 V N_f) * sum_{p} log(2 - 4 lambda_N - 4*kappa*sum{cos(P_mu)} + 8 lambda_N*(m^2+s^2) + 18*lambda_6_N(m^4 + s^4 + 6 m^2 s^2))  excluding zero and staggered mode
	//using the shortcut: D_p=2 - 4 lambda_N - 4*kappa*sum{cos(P_mu)} + 8 lambda_N*(m^2+s^2) + 18*lambda_6_N(m^4 + s^4 + 6 m^2 s^2)
	//with that:
	//dU/dm += -1/(2 V N_f) * sum_{p}[ (16*lambda_N*m + 72*lambda_6_N*(m^3 + 3*m*s^2)) / D_p ]
	//dU/ds += -1/(2 V N_f) * sum_{p}[ (16*lambda_N*s + 72*lambda_6_N*(s^3 + 3*m^2*s)) / D_p ]
	//
	//NOTE, not both possibilities can be used at the same time
	if(useBosonicLoop && useImprovedGaussian)
	{
		std::cerr <<"Error, using the bosonic loop and the improved gaussian contribution at the same time is not implemented" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(useBosonicLoop && lambda_6_N!=0.0)
	{
		std::cerr <<"Error, using the bosonic loop non-zero lambda_6_N at the same time is not implemented" <<std::endl;
		exit(EXIT_FAILURE);
	}
	dU_ov_dm=0.0; dU_ov_ds=0.0;
	computeFermionicContribution_onlyGradient_FromStoredEigenvalues(magnetization, staggeredMagnetization, dU_ov_dm, dU_ov_ds);
	if(useBosonicLoop && !bosonicLoopSet){ bosonicLoop=computeBosonicPropagatorSum_fromStoredSumOfCos();}
	double mSq(magnetization*magnetization), sSq(staggeredMagnetization*staggeredMagnetization);
	//std::cout <<"ferm. contr to grad: dU_f/dm=" <<dU_ov_dm <<"   dU_f/ds=" <<dU_ov_ds <<std::endl;
	dU_ov_dm += -16.0*kappa_N*magnetization + 2.0*magnetization;
	dU_ov_dm += lambda_N*( 4.0*mSq*magnetization + 12.0*magnetization*sSq - 4.0*magnetization);
	dU_ov_dm += lambda_6_N*( 6.0*mSq*mSq*magnetization + 60.0*mSq*magnetization*sSq + 30.0*magnetization*sSq*sSq);
	if(useBosonicLoop){ dU_ov_dm += 32.0*lambda_N*magnetization*bosonicLoop/static_cast< double >(N_f); }
	
	dU_ov_ds += +16.0*kappa_N*staggeredMagnetization + 2.0*staggeredMagnetization;
	dU_ov_ds += lambda_N*( 4.0*sSq*staggeredMagnetization + 12.0*mSq*staggeredMagnetization - 4.0*staggeredMagnetization);
	dU_ov_ds += lambda_6_N*( 6.0*sSq*sSq*staggeredMagnetization + 30.0*mSq*mSq*staggeredMagnetization + 60.0*mSq*sSq*staggeredMagnetization);
	if(useBosonicLoop){ dU_ov_ds += 32.0*lambda_N*staggeredMagnetization*bosonicLoop/static_cast< double >(N_f); }
	
	if(useImprovedGaussian)
	{
		double bosDet_ov_dm(0.0), bosDet_ov_ds(0.0);
		computeBosonicDeterminantContributionForImprovedGaussian_onlyGradient_fromStoredSumOfCos(magnetization, staggeredMagnetization, bosDet_ov_dm, bosDet_ov_ds);
		dU_ov_dm+=bosDet_ov_dm;
		dU_ov_ds+=bosDet_ov_ds;
	}
	//std::cout <<"full grad: dU_f/dm=" <<dU_ov_dm <<"   dU_f/ds=" <<dU_ov_ds <<std::endl;

}




void constrainedEffectivePotential::fermionicContributionInline_onlyGradient_FromStoredEigenvalues(int index, double ySq_mSq, double ySq_sSq, double &dUf_ov_dm, double &dUf_ov_ds)
{
	double a,b,one_ov_A;
	a = absNuP[index]*absNuVarP[index] + (ySq_mSq - ySq_sSq)*absGammaP[index]*absGammaVarP[index];
	b = absGammaP[index]*absNuVarP[index] - absNuP[index]*absGammaVarP[index];
	one_ov_A=1.0/(a*a + ySq_mSq*b*b);
	//std::cout <<"one_ov_A =" <<one_ov_A <<std::endl;
	dUf_ov_dm=(2.0*absGammaP[index]*absGammaVarP[index]*a + b*b)*one_ov_A*factorOfMomentum[index];
	dUf_ov_ds=(absGammaP[index]*absGammaVarP[index]*a)*one_ov_A*factorOfMomentum[index];
}







void constrainedEffectivePotential::computeFermionicContribution_onlyGradient_FromStoredEigenvalues( const double magnetization, const double staggeredMagnetization, double &dUf_ov_dm, double &dUf_ov_ds)
{
	//U_f = -2/V*sum log[ a^2 + y^2*m^2*b^2 ] = -2/v*sum log[A]
	//a = a(m,s) = |nu(p)|*|nu(varP)| + y^2*(m^2-s^2)*|gamma(p)|*|gamma(varP)|
	//b = |gamma(p)|*|nu(varP)| - |nu(p)|*|gamma(varP)|
	//A = a^2 + y^2*m^2*b^2
	//dU_f/dm = -2/V*sum[ (4*y^2*m * |gamma(p)|*|gamma(varP)|*a + 2*m*y^2*b^2)/(a^2 + y^2*m^2*b^2 ) ] = -4*y^2*m/V * sum[ (2* |gamma(p)|*|gamma(varP)|*a + b^2)/(a^2 + y^2*m^2*b^2 ) ]
	//dU_f/dm = -4*y^2*m/V * sum[ (2* |gamma(p)|*|gamma(varP)|*a + b^2)/A ]
	//dU_f/ds = -2/V*sum[ (-4*y^2*s * |gamma(p)|*|gamma(varP)|*a )/(a^2 + y^2*m^2*b^2 ) ] = +8*y^2*s/V * sum[ ( |gamma(p)|*|gamma(varP)|*a )/(a^2 + y^2*m^2*b^2 ) ]
	//dU_f/ds = 8*y^2*m/V * sum[ ( |gamma(p)|*|gamma(varP)*a|)/A ]
	if(yukawa_N < 0.0)
	{
		std::cerr <<"Error, no yukawa coupling set in constrainedEffectivePotential::computeFermionicContribution_onlyGradient_FromStoredEigenvalues( const double magnetization, const double staggeredMagnetization, double &dUf_ov_dm, double &dUf_ov_ds)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	dUf_ov_dm=0.0; dUf_ov_ds=0.0; 
	if(yukawa_N==0.0){ return; }
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
	double ySquared(yukawa_N*yukawa_N);
	double ySquaredmSquared(ySquared*mSq);
	double ySquaredsSquared(ySquared*sSq);
	
	double dummyForAddition_dm(0.0), dummyForAddition_ds(0.0);
	
	for(int index=0; index<numberOfDistingtMomenta; ++index)
	{
		fermionicContributionInline_onlyGradient_FromStoredEigenvalues(index, ySquaredmSquared, ySquaredsSquared, dummyForAddition_dm, dummyForAddition_ds);
		dUf_ov_dm+=dummyForAddition_dm;
		dUf_ov_ds+=dummyForAddition_ds;
	}
	
	dUf_ov_dm*=-4.0*ySquared*magnetization;
	dUf_ov_ds*=8.0*ySquared*staggeredMagnetization;
	
	dUf_ov_dm/=static_cast< double >(L0); dUf_ov_dm/=static_cast< double >(L1); dUf_ov_dm/=static_cast< double >(L2); dUf_ov_dm/=static_cast< double >(L3);
	dUf_ov_ds/=static_cast< double >(L0); dUf_ov_ds/=static_cast< double >(L1); dUf_ov_ds/=static_cast< double >(L2); dUf_ov_ds/=static_cast< double >(L3);
	
	return;
}



void constrainedEffectivePotential::computeBosonicDeterminantContributionForImprovedGaussian_onlyGradient_fromStoredSumOfCos(const double magnetization, const double staggeredMagnetization, double &dU_ov_dm, double &dU_ov_ds)
{
	//computes:
	//using the shortcut: D_p=2 - 4 lambda_N - 4*kappa*sum{cos(P_mu)} + 8 lambda_N*(m^2+s^2) + 18*lambda_6_N(m^4 + s^4 + 6 m^2 s^2)
	//with that:
	//dU_bosDet/dm += -1/(2 V N_f) * sum_{p}[ (16*lambda_N*m + 72*lambda_6_N*(m^3 + 3*m*s^2)) / D_p ]
	//dU_bosDet/ds += -1/(2 V N_f) * sum_{p}[ (16*lambda_N*s + 72*lambda_6_N*(s^3 + 3*m^2*s)) / D_p ]
	double mSq(magnetization*magnetization), sSq(staggeredMagnetization*staggeredMagnetization);
	double constantPart(2.0 - 4.0*lambda_N + 8.0*lambda_N*( mSq + sSq ) + 18.0*lambda_6_N*(mSq*mSq + sSq*sSq + 6.0*mSq*sSq) );
	double fourKappa=4.0*kappa_N;
	
	dU_ov_dm=0.0;
	dU_ov_ds=0.0;
	double dummyForAddition(0.0);
	
	for(int index=0; index<numberOfDistingtMomenta_bosonic; ++index)
	{
		dummyForAddition += factorOfMomentum_bosonic[index]/(constantPart - fourKappa*sumOfCosOfPmu[index]);
	}
	
	dummyForAddition*=-0.5;
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
	dummyForAddition/=static_cast< double >(N_f);
	
	dU_ov_dm=dummyForAddition*(16.0*lambda_N*magnetization + 72.0*lambda_6_N*(mSq*magnetization + 3.0*magnetization*sSq));
	dU_ov_ds=dummyForAddition*(16.0*lambda_N*staggeredMagnetization + 72.0*lambda_6_N*(sSq*staggeredMagnetization + 3.0*mSq*staggeredMagnetization));
	
}


void constrainedEffectivePotential::computeConstrainedEffectivePotential_onlyGradient_gsl(const gsl_vector *mags, gsl_vector *gradient_of_U)
{
	double magnetization=gsl_vector_get(mags,0);
	double staggeredMagnetization=gsl_vector_get(mags,1);
	
	double dU_dm(0.0), dU_ds(0.0);
	
	computeConstrainedEffectivePotential_onlyGradient( magnetization, staggeredMagnetization, dU_dm, dU_ds);
	gsl_vector_set(gradient_of_U, 0, dU_dm);
	gsl_vector_set(gradient_of_U, 1, dU_ds);
	
}


void wrapper_computeConstrainedEffectivePotential_onlyGradient_gsl(const gsl_vector *mags, void *params, gsl_vector *gradient_of_U)
{
	constrainedEffectivePotential *CEP = (constrainedEffectivePotential *)params;
	CEP->computeConstrainedEffectivePotential_onlyGradient_gsl(mags, gradient_of_U);
}


/*
void constrainedEffectivePotential::computeFermionicContribution_onlyGradient_qad( const double magnetization, const double staggeredMagnetization, double &dUf_ov_dm, double &dUf_ov_ds)
{
	//U_f = -2/V*sum log[ a^2 + y^2*m^2*b^2 ] = -2/v*sum log[A]
	//a = a(m,s) = |nu(p)|*|nu(varP)| + y^2*(m^2-s^2)*|gamma(p)|*|gamma(varP)|
	//b = |gamma(p)|*|nu(varP)| - |nu(p)|*|gamma(varP)|
	//A = a^2 + y^2*m^2*b^2
	//dU_f/dm = -2/V*sum[ (4*y^2*m * |gamma(p)|*|gamma(varP)|*a + 2*m*y^2*b^2)/(a^2 + y^2*m^2*b^2 ) ] = -4*y^2*m/V * sum[ (2* |gamma(p)|*|gamma(varP)|*a + b^2)/(a^2 + y^2*m^2*b^2 ) ]
	//dU_f/dm = -4*y^2*m/V * sum[ (2* |gamma(p)|*|gamma(varP)|*a + b^2)/A ]
	//dU_f/ds = -2/V*sum[ (-4*y^2*s * |gamma(p)|*|gamma(varP)|*a )/(a^2 + y^2*m^2*b^2 ) ] = +8*y^2*s/V * sum[ ( |gamma(p)|*|gamma(varP)|*a )/(a^2 + y^2*m^2*b^2 ) ]
	//dU_f/ds = 8*y^2*m/V * sum[ ( |gamma(p)|*|gamma(varP)*a|)/A ]
	
	if(yukawa_N < 0.0)
	{
		std::cerr <<"Error, no yukawa coupling set in constrainedEffectivePotential::computeFermionicContribution_onlyGradient_qad(double magSquared, double stagMagSquared)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	dUf_ov_dm=0.0; dUf_ov_ds=0.0; 
	if(yukawa_N==0.0){ return; }
	
// 	double one_ov_twoRho=0.5/rho;
	
	double ySquared(yukawa_N*yukawa_N);
	double ySq_times_mSq_mi_sSq( ySquared*(magnetization*magnetization - staggeredMagnetization*staggeredMagnetization) );
	double ySq_times_mSq( ySquared*magnetization*magnetization );
	
	std::complex< double > nuOfP(0.,0.), nuOfVarP(0.,0.);
	double abs_nuOfP, abs_nuOfVarP;
	double abs_gammaOfP, abs_gammaOfVarP;
	double a,b,one_ov_A;
	
	for(int l0=0; l0<L0; ++l0)
	{
		for(int l1=0; l1<L1; ++l1)
		{
			for(int l2=0; l2<L2; ++l2)
			{
				for(int l3=0; l3<L3; ++l3)
				{
					//NOTE since only large small lattices will be used, at some point, make an array of abs(nu)
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					
					a = abs_nuOfP*abs_nuOfVarP + ySq_times_mSq_mi_sSq*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					one_ov_A=1.0/(a*a + ySq_times_mSq*b*b);
					dUf_ov_dm+=(2.0*abs_gammaOfP*abs_gammaOfVarP*a + b*b)*one_ov_A;
					dUf_ov_ds+=(abs_gammaOfP*abs_gammaOfVarP*a)*one_ov_A;
					
				}
			}
		}
	}
	dUf_ov_dm*=-4.0*ySquared*magnetization;
	dUf_ov_ds*=8.0*ySquared*staggeredMagnetization;
	
	dUf_ov_dm/=static_cast< double >(L0); dUf_ov_dm/=static_cast< double >(L1); dUf_ov_dm/=static_cast< double >(L2); dUf_ov_dm/=static_cast< double >(L3);
	dUf_ov_ds/=static_cast< double >(L0); dUf_ov_ds/=static_cast< double >(L1); dUf_ov_ds/=static_cast< double >(L2); dUf_ov_ds/=static_cast< double >(L3);
}
*/