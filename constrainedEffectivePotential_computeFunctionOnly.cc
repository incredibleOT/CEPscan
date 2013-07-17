#include "constrainedEffectivePotential.h"


double constrainedEffectivePotential::computeConstrainedEffectivePotential_onlyFunction(const double magnetization, const double staggeredMagnetization)
{
	//eq 4.35 philipp's thesis
	//U= -8*kappa*(m^2-s^2) + m^2 + s^2 + lambda*( m^4 + s^4 + 6*m^2s^2 - 2*(m^2+s^2)) + lambda_6(m^6 + s^6 + 15*m^4*s^2 + 15*m^2*s^4) U_f(m,s)
	//dU/dm = dU_f/d_m - 16 kappa * m + 2*m + lambda*(4*m^3 + 12 *m*s^2 - 4*m) + lambda_6*(6*m^5 + 60*m^3*s^2 + 30*m*s^4)
	//dU/ds = dU_f/d_m + 16 kappa * s + 2*s + lambda*(4*s^3 + 12 *m^2*s - 4*s) + lambda_6*(6*s^5 + 30*m^4*s + 60*m^2*s^3)
	//
	//NOTE there are three possible improvements for the potential
	//
	//the 1st oneis a term coming from the boson loop, or better: the first order in lambda:
	//U+=16*(m^2+s^2)*lambda_N*N_f^{-1}*P_B_naive
	//P_B_naive= 1/V * sum_{p} 1/(2-4*lambda_N-4*kappa*sum{cos(P_mu)}) excluding zero and staggered mode
	//NOTE, does not work with lambda_6!=0
	//
	//the 2nd possibility is to include an improved bosonic determinant:
	//U+=-1/(2 V N_f) * sum_{p} log(2 - 4 lambda_N - 4*kappa*sum{cos(P_mu)} + 8 lambda_N*(m^2+s^2) + 18*lambda_6_N(m^4 + s^4 + 6 m^2 s^2))  excluding zero and staggered mode
	//in short: U+=-1/(2 V N_f)*sum_{p} D_p(m,s)
	//with: D_p(m,s)=2 - 4 lambda_N - 4*kappa*sum{cos(P_mu)} + 8 lambda_N*(m^2+s^2) + 18*lambda_6_N(m^4 + s^4 + 6 m^2 s^2)
	//NOTE, not both of the possibilities above can be used at the same time
	//
	//3rd: Also consider the first order in lambda and lambda_6 (makes only sence, if the improved tree-level is used)
	//U+=8/N_f^2*(lambda_N + lambda_6N * (9*(m^2+s^2)))* P_B^2 + 24*lambda_6N/N_f^3*P_B^3
	//with P_B=1/V*sum_{p} 1/D_p(m,s)   excluding the zero and staggered mode
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
	
	double result(0.0);
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
	if(useBosonicLoop && !bosonicLoopSet){ bosonicLoop=computeBosonicPropagatorSum_fromStoredSumOfCos();}
	
	double fermionicContribution_sumOfLogs=computeFermionicContribution_onlyFunction_FromStoredEigenvalues(magnetization,staggeredMagnetization);
	result+=fermionicContribution_sumOfLogs;
	result+= mSq + sSq - 8.0 * kappa_N * (mSq - sSq);
	result+= lambda_N * ( mSq*mSq + sSq*sSq + 6.0*mSq*sSq -2.0*(mSq + sSq) );
	result+= lambda_6_N*( mSq*mSq*mSq + sSq*sSq*sSq + 15.0*mSq*mSq*sSq + 15.0*mSq*sSq*sSq );
	if(useBosonicLoop){ result+= 16.0*lambda_N*( mSq + sSq )*bosonicLoop/static_cast< double >(N_f); }
	
	if(useImprovedGaussian && !useImprovedFirstOrder)
	{
		double bosDet( computeBosonicDeterminantContributionForImprovedGaussian_onlyFunction_fromStoredSumOfCos(magnetization, staggeredMagnetization) );
		result+=bosDet;
	}
	else if(useImprovedFirstOrder)
	{
		double bosDetAnd1stOrder( computeImprovedBosDetAndFirstOrderContribution_onlyFunction_fromStoredSumOfCos(magnetization, staggeredMagnetization) );
		result+=bosDetAnd1stOrder;
	}
	if(result!=result){result=-log(0.0);}
	return result;
}





double constrainedEffectivePotential::fermionicContributionInline_onlyFunction_FromStoredEigenvalues(int index, double ySq_mSq, double ySq_sSq)
{
	
	double a,b;
	a = absNuP[index]*absNuVarP[index] + (ySq_mSq - ySq_sSq)*absGammaP[index]*absGammaVarP[index];
	b = absGammaP[index]*absNuVarP[index] - absNuP[index]*absGammaVarP[index];
	return factorOfMomentum[index]*log(a*a + ySq_mSq*b*b);
}



double constrainedEffectivePotential::computeFermionicContribution_onlyFunction_FromStoredEigenvalues(const double magnetization, const double staggeredMagnetization)
{
	if(yukawa_N < 0.0)
	{
		std::cerr <<"Error, no yukawa coupling set in constrainedEffectivePotential::computeFermionicContribution_onlyFunction_qad(double magnetization, double staggeredMagnetization)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(yukawa_N==0.0){ return 0.0; }
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
	double ySquared(yukawa_N*yukawa_N);
	double ySquaredmSquared(ySquared*mSq);
	double ySquaredsSquared(ySquared*sSq);
	double dummyForAddition(0.0);
	
	for(int index=0; index<numberOfDistingtMomenta; ++index)
	{
		dummyForAddition+=fermionicContributionInline_onlyFunction_FromStoredEigenvalues(index, ySquaredmSquared, ySquaredsSquared);
	}
	
	dummyForAddition*=-2.0;
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
	return dummyForAddition;
}


double constrainedEffectivePotential::computeConstrainedEffectivePotential_onlyFunction_gsl(const gsl_vector *mags)
{
	double magnetization=gsl_vector_get(mags,0);
	double staggeredMagnetization=gsl_vector_get(mags,1);
	
	return computeConstrainedEffectivePotential_onlyFunction(magnetization, staggeredMagnetization);
}



double constrainedEffectivePotential::computeBosonicDeterminantContributionForImprovedGaussian_onlyFunction_qad(const double magnetization, const double staggeredMagnetization)
{
	//computes: -1/(2 V N_f) * sum_{p} log(2 - 4 lambda_N - 4*kappa*sum{cos(P_mu)} + 8 lambda_N*(m^2+s^2) + 18*lambda_6_N(m^4 + s^4 + 6 m^2 s^2))  excluding zero and staggered mode
	double result(0.0);
	double mSq(magnetization*magnetization), sSq(staggeredMagnetization*staggeredMagnetization);
	double constantPart(2.0 - 4.0*lambda_N + 8.0*lambda_N*( mSq + sSq ) + 18.0*lambda_6_N*(mSq*mSq + sSq*sSq + 6.0*mSq*sSq) );
	double fourKappa=4.0*kappa_N;
	
	int L0_half=L0/2, L1_half=L1/2, L2_half=L2/2, L3_half=L3/2;
	for(int l0=0; l0<L0; ++l0)
	{
		for(int l1=0; l1<L1; ++l1)
		{
			for(int l2=0; l2<L2; ++l2)
			{
				for(int l3=((l0+l1+l2)?0:1); l3<L3_half; ++l3)
				{
					
					result+=log(constantPart - fourKappa*(cosOfPmu_L0[l0] + cosOfPmu_L1[l1] + cosOfPmu_L2[l2] + cosOfPmu_L3[l3]));
// 					counter++;
				}
				for(int l3=((l0==L0_half && l1==L1_half && l2==L2_half )?L3_half+1:L3_half); l3<L3; ++l3)
				{
					result+=log(constantPart - fourKappa*(cosOfPmu_L0[l0] + cosOfPmu_L1[l1] + cosOfPmu_L2[l2] + cosOfPmu_L3[l3]));
// 					counter++;
				}
				
			}
		}
	}
	result*=-0.5;
	result/=static_cast< double >(L0); result/=static_cast< double >(L1); result/=static_cast< double >(L2); result/=static_cast< double >(L3); 
	result/=static_cast< double >(N_f);
	return result;
}




double constrainedEffectivePotential::computeBosonicDeterminantContributionForImprovedGaussian_onlyFunction_fromStoredSumOfCos(const double magnetization, const double staggeredMagnetization)
{
	//computes: -1/(2 V N_f) * sum_{p} log(2 - 4 lambda_N - 4*kappa*sum{cos(P_mu)} + 8 lambda_N*(m^2+s^2) + 18*lambda_6_N(m^4 + s^4 + 6 m^2 s^2))  excluding zero and staggered mode
	double result(0.0);
	double mSq(magnetization*magnetization), sSq(staggeredMagnetization*staggeredMagnetization);
	double constantPart(2.0 - 4.0*lambda_N + 8.0*lambda_N*( mSq + sSq ) + 18.0*lambda_6_N*(mSq*mSq + sSq*sSq + 6.0*mSq*sSq) );
	double fourKappa=4.0*kappa_N;
	
	for(int index=0; index<numberOfDistingtMomenta_bosonic; ++index)
	{
		result+= factorOfMomentum_bosonic[index]*log(constantPart - fourKappa*sumOfCosOfPmu[index]);
	}
	
	result*=-0.5;
	result/=static_cast< double >(L0); result/=static_cast< double >(L1); result/=static_cast< double >(L2); result/=static_cast< double >(L3); 
	result/=static_cast< double >(N_f);
	return result;
}



double constrainedEffectivePotential::computeImprovedBosDetAndFirstOrderContribution_onlyFunction_fromStoredSumOfCos(const double magnetization, const double staggeredMagnetization)
{
	double U_BosDet(0.0),U_1st(0.0);
	computeImprovedBosDetAndFirstOrderContribution_onlyFunction_fromStoredSumOfCos(magnetization, staggeredMagnetization, U_BosDet, U_1st);
	return U_BosDet+U_1st;
}


void constrainedEffectivePotential::computeImprovedBosDetAndFirstOrderContribution_onlyFunction_fromStoredSumOfCos(const double magnetization, const double staggeredMagnetization, double &U_BosDet, double &U_1st)
{
	//computes:
	//contribution from the bosonic determinant in case of the improved tree-level and the first order contribution
	//
	//U_BosDet = - 1/(2 V N_f) * sum_{p} log(D_p(m,s) 
	//and
	//U_1st =  8/N_f^2*(lambda_N + lambda_6N * (9*(m^2+s^2)))* P_B^2 + 24*lambda_6N/N_f^3*P_B^3
	//with 
	//P_B=\sum_{p} D_p(m,s)   excluding the zero and staggered mode
	//and D_p(m,s)=2 - 4 lambda_N - 4*kappa*sum{cos(P_mu)} + 8 lambda_N*(m^2+s^2) + 18*lambda_6_N(m^4 + s^4 + 6 m^2 s^2)
	U_BosDet=0.0;
	U_1st=0.0;
	
	double mSq(magnetization*magnetization), sSq(staggeredMagnetization*staggeredMagnetization);
	double constantPart(2.0 - 4.0*lambda_N + 8.0*lambda_N*( mSq + sSq ) + 18.0*lambda_6_N*(mSq*mSq + sSq*sSq + 6.0*mSq*sSq) );
	double fourKappa=4.0*kappa_N;
	
	double dummyForLog(0.0), dummyForAddition(0.0);
	for(int index=0; index<numberOfDistingtMomenta_bosonic; ++index)
	{
		dummyForLog+= factorOfMomentum_bosonic[index]*log(constantPart - fourKappa*sumOfCosOfPmu[index]);
		dummyForAddition+=factorOfMomentum_bosonic[index]/(constantPart - fourKappa*sumOfCosOfPmu[index]);
	}
	
	dummyForLog*=-0.5;
	dummyForLog/=static_cast< double >(L0); dummyForLog/=static_cast< double >(L1); dummyForLog/=static_cast< double >(L2); dummyForLog/=static_cast< double >(L3); 
	dummyForLog/=static_cast< double >(N_f);
	
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
	U_BosDet=dummyForLog;
	U_1st=8.0/static_cast< double >(N_f*N_f)*(lambda_N + lambda_6_N*9.0*(mSq+sSq))*dummyForAddition*dummyForAddition;
	U_1st+=24.0*lambda_6_N/static_cast< double >(N_f*N_f*N_f)*dummyForAddition*dummyForAddition*dummyForAddition;
}



double wrapper_computeConstrainedEffectivePotential_onlyFunction_gsl(const gsl_vector *mags, void *params)
{
	constrainedEffectivePotential *CEP = (constrainedEffectivePotential *)params;
	return CEP->computeConstrainedEffectivePotential_onlyFunction_gsl(mags);
}




/*
//computation of the fermionic contribution
//just computes the sum of logs as in philipp's thesis 4.35
double constrainedEffectivePotential::computeFermionicContribution_onlyFunction_qad(const double magnetization, const double staggeredMagnetization)
{
// computes U_f
// U_f=-1/V * sum_p log( A^2 ) = -2/V *sum_p log( A ) 
// A = a^2 + m^2*y^2* b^2
// a = |nu(p)|*|nu(varP)| + y_N^2*(m^2-s^2)*|gamma(p)|*|gamma(varP)|
// b = |gamma(p)|*|nu(varP)| - |nu(p)|*|gamma(varP)|
// gamma(p) = 1 - nu(p)/2*rho
// 
	if(yukawa_N < 0.0)
	{
		std::cerr <<"Error, no yukawa coupling set in constrainedEffectivePotential::computeFermionicContribution_onlyFunction_qad(double magnetization, double staggeredMagnetization)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(yukawa_N==0.0){ return 0.0; }
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
// 	const double two_PI(atan(1) * 8.0);
// 	double one_ov_L0(1.0/L0);
// 	double one_ov_L1(1.0/L1);
// 	double one_ov_L2(1.0/L2);
// 	double one_ov_L3(1.0/L3);
// 	double toAddForL3 = antiperiodicBC_L3 ? two_PI/(2.0*L3) : 0.0;
// 	double one_ov_twoRho=0.5/rho;
	double ySquared(yukawa_N*yukawa_N);
	
	//dummies needed for the summation
	double dummyForAddition(0.0);
	std::complex< double > nuOfP(0.,0.), nuOfVarP(0.,0.);
	double abs_nuOfP, abs_nuOfVarP;
	double abs_gammaOfP, abs_gammaOfVarP;
	double a,b;
	for(int l0=0; l0<L0; ++l0)
	{
// 		double p0 = two_PI * l0 * one_ov_L0 ; //2*pi*n/L
// 		double varp0 = two_PI * ((l0+L0/2)%L0) * one_ov_L0;
		for(int l1=0; l1<L1; ++l1)
		{
// 			double p1 = two_PI * l1 * one_ov_L1 ;
// 			double varp1 = two_PI * ((l1+L1/2)%L1) * one_ov_L1;
			for(int l2=0; l2<L2; ++l2)
			{
// 				double p2 = two_PI * l2 * one_ov_L2 ;
// 				double varp2 = two_PI * ((l2+L2/2)%L2) * one_ov_L2;
				for(int l3=0; l3<L3; ++l3)
				{
// 					double p3 = two_PI * l3 * one_ov_L3 + toAddForL3;
// 					double varp3 = two_PI * ((l3+L3/2)%L3) * one_ov_L3 + toAddForL3;
					//NOTE since only large small lattices will be used, at some point, make an array of abs(nu)
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
// 					std::cout <<"l0=" <<l0 <<"  l1=" <<l1 <<"  l2=" <<l2 <<" l3=" <<l3; 
// 					std::cout <<":   nu(p)=" <<nuOfP <<"   abs(nu(p))=" << abs_nuOfP
// 					<<"   nu(var(p))="<<nuOfVarP <<"   abs(nu(varP))=" << abs_nuOfVarP <<std::endl;
// 					std::cout <<"   abs(nu(p))=" << abs_nuOfP <<"   abs(nu(varP))=" << abs_nuOfVarP 
// 					          <<"   abs(gamma(p)=" <<abs_gammaOfP <<"   abs(gamma(varP)=" <<abs_gammaOfVarP << std::endl;
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=log(a*a + ySquared*mSq*b*b);
					//std::cout <<"dummyForAddition=" <<dummyForAddition <<std::endl;
				}
			}
		}
	}
	dummyForAddition*=-2.0;
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
	return dummyForAddition;
}
*/



/*
double constrainedEffectivePotential::computeFermionicContribution_onlyFunction_LcubeTimesLt(const double magnetization, const double staggeredMagnetization)
{
// computes U_f
// U_f=-1/V * sum_p log( A^2 ) = -2/V *sum_p log( A ) 
// A = a^2 + m^2*y^2* b^2
// a = |nu(p)|*|nu(varP)| + y_N^2*(m^2-s^2)*|gamma(p)|*|gamma(varP)|
// b = |gamma(p)|*|nu(varP)| - |nu(p)|*|gamma(varP)|
// gamma(p) = 1 - nu(p)/2*rho
//
//Function is speciallized for an L^3xL_t lattice (uses permutation of components)
//also abuses, that the contributions of p_mu=pi-x is the same as for p_mu=pi+x
//cases that can occur (notation: p,q,r \in [1,L/2-1]:,p<q<r) 
//timecomponent is not shown, gives an additional factor of 2 for values not equal zero or pi (for even BC)
//for antiperiodic BC in time, every contribution gives a factor 2
// momentum    factor
//
// p,q,r         48
// p,q,L/2       24
// 0,p,q         24
// 0,p,L/2       12
//
// p,p,q         24
// p,q,q         24
// p,p,L/2       12
// p,L/2,L/2     6
// 0,0,p         6
// 0,p,p         12
// 0,L/2,L/2     3
// 0,0,L/2       3
//
// p,p,p         8
// L/2,L/2,L/2   1
// 0,0,0         1




	if(yukawa_N < 0.0)
	{
		std::cerr <<"Error, no yukawa coupling set in constrainedEffectivePotential::computeFermionicContribution_onlyFunction_qad(double magnetization, double staggeredMagnetization)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(yukawa_N==0.0){ return 0.0; }
	if(!(L0==L1 && L0==L2))
	{
		std::cerr <<"Error, no L^3xL_t lattice given in constrainedEffectivePotential::computeFermionicContribution_onlyFunction_LcubeTimesLt" <<std::endl;
		exit(EXIT_FAILURE);
	}
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
// 	double one_ov_twoRho=0.5/rho;
	double ySquared(yukawa_N*yukawa_N);
	
	//dummies needed for the summation
	double dummyForAddition(0.0);
	std::complex< double > nuOfP(0.,0.), nuOfVarP(0.,0.);
	double abs_nuOfP, abs_nuOfVarP;
	double abs_gammaOfP, abs_gammaOfVarP;
	double a,b;
	int Lhalf=L0/2;
	int LT_half=L3/2;
	for(int l3=1; l3<LT_half; ++l3)
	{
		for(int l0=1; l0<Lhalf; ++l0)
		{
			for(int l1=l0+1; l1<Lhalf; ++l1)
			{
				for(int l2=l1+1; l2<Lhalf; ++l2)
				{
					//p,q,r 48
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=48.0*log(a*a + ySquared*mSq*b*b);
				}
				{
					int l2=l1;
					//p,q,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=24.0*log(a*a + ySquared*mSq*b*b);
					l2=l0;
					//p,p,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=24.0*log(a*a + ySquared*mSq*b*b);
					l2=Lhalf;
					//p,q,L/2 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=24.0*log(a*a + ySquared*mSq*b*b);
					l2=0;
					//0,p,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=24.0*log(a*a + ySquared*mSq*b*b);
				}
			}
			{
				int l1=l0;
				{
					int l2=l1;
					//p,p,p, 8
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=8.0*log(a*a + ySquared*mSq*b*b);
					l2=Lhalf;
					//p,p,L/2 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=12.0*log(a*a + ySquared*mSq*b*b);
					l2=0;
					//0,p,p, 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=12.0*log(a*a + ySquared*mSq*b*b);
				}
				l1=Lhalf;
				{
					int l2=Lhalf;
					//p,L/2,L/2, 6
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=6.0*log(a*a + ySquared*mSq*b*b);
					l2=0;
					//0,p,L/2, 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=12.0*log(a*a + ySquared*mSq*b*b);
				}
				l1=0;
				{
					int l2=0;
					//0,0,p  6
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=6.0*log(a*a + ySquared*mSq*b*b);
				}
			}
		}
		{
			int l0=Lhalf;
			{
				int l1=Lhalf;
				{
					int l2=Lhalf;
					//L/2,L/2,L/2, 1
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=1.0*log(a*a + ySquared*mSq*b*b);
					l2=0;
					//0,L/2,L/2, 3
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=3.0*log(a*a + ySquared*mSq*b*b);
				}
				l1=0;
				{
					int l2=0;
					//0,0,L/2, 3
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=3.0*log(a*a + ySquared*mSq*b*b);
				}
			}
			l0=0;
			{
				int l1=0;
				{
					int l2=0;
					//0,0,0, 1
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=1.0*log(a*a + ySquared*mSq*b*b);
				}
			}
		}
	}
	//now multiplzy with two, due to the symmetry in L3
	dummyForAddition*=2.0;
	//now do the same thing for l3=0 (and l3=Lt_half if periodic BC are used)
	double BC_factor=(antiperiodicBC_L3 ? 2.0 :1.0);
	for(int l3=0; l3<=LT_half; l3+=(antiperiodicBC_L3 ? LT_half*2 :LT_half))
	{
		for(int l0=1; l0<Lhalf; ++l0)
		{
			for(int l1=l0+1; l1<Lhalf; ++l1)
			{
				for(int l2=l1+1; l2<Lhalf; ++l2)
				{
					//p,q,r 48
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=48.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
				}
				{
					int l2=l1;
					//p,q,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=24.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
					l2=l0;
					//p,p,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=24.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
					l2=Lhalf;
					//p,q,L/2 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=24.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
					l2=0;
					//0,p,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=24.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
				}
			}
			{
				int l1=l0;
				{
					int l2=l1;
					//p,p,p, 8
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=8.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
					l2=Lhalf;
					//p,p,L/2 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=12.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
					l2=0;
					//0,p,p, 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=12.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
				}
				l1=Lhalf;
				{
					int l2=Lhalf;
					//p,L/2,L/2, 6
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=6.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
					l2=0;
					//0,p,L/2, 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=12.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
				}
				l1=0;
				{
					int l2=0;
					//0,0,p  6
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=6.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
				}
			}
		}
		{
			int l0=Lhalf;
			{
				int l1=Lhalf;
				{
					int l2=Lhalf;
					//L/2,L/2,L/2, 1
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=1.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
					l2=0;
					//0,L/2,L/2, 3
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=3.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
				}
				l1=0;
				{
					int l2=0;
					//0,0,L/2, 3
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=3.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
				}
			}
			l0=0;
			{
				int l1=0;
				{
					int l2=0;
					//0,0,0, 1
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					abs_nuOfP=std::abs(nuOfP); 
					abs_nuOfVarP=std::abs(nuOfVarP);
					abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=1.0*log(a*a + ySquared*mSq*b*b)*BC_factor;
				}
			}
		}
	}
	//now multiply with two, if antiperiodic BC in time are chosen, due to the symmetry in L3
// 	if(antiperiodicBC_L3){ dummyForAddition*=2.0; }//NOTE, this is wrong!
	dummyForAddition*=-2.0;
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
	return dummyForAddition;
}




double constrainedEffectivePotential::computeFermionicContribution_onlyFunction_LcubeTimesLt_withInline(const double magnetization, const double staggeredMagnetization)
{
// computes U_f
// U_f=-1/V * sum_p log( A^2 ) = -2/V *sum_p log( A ) 
// A = a^2 + m^2*y^2* b^2
// a = |nu(p)|*|nu(varP)| + y_N^2*(m^2-s^2)*|gamma(p)|*|gamma(varP)|
// b = |gamma(p)|*|nu(varP)| - |nu(p)|*|gamma(varP)|
// gamma(p) = 1 - nu(p)/2*rho
//
//Function is speciallized for an L^3xL_t lattice (uses permutation of components)
//also abuses, that the contributions of p_mu=pi-x is the same as for p_mu=pi+x
//cases that can occur (notation: p,q,r \in [1,L/2-1]:,p<q<r) 
//timecomponent is not shown, gives an additional factor of 2 for values not equal zero or pi (for even BC)
//for antiperiodic BC in time, every contribution gives a factor 2
// momentum    factor
//
// p,q,r         48
// p,q,L/2       24
// 0,p,q         24
// 0,p,L/2       12
//
// p,p,q         24
// p,q,q         24
// p,p,L/2       12
// p,L/2,L/2     6
// 0,0,p         6
// 0,p,p         12
// 0,L/2,L/2     3
// 0,0,L/2       3
//
// p,p,p         8
// L/2,L/2,L/2   1
// 0,0,0         1




	if(yukawa_N < 0.0)
	{
		std::cerr <<"Error, no yukawa coupling set in constrainedEffectivePotential::computeFermionicContribution_onlyFunction_qad(double magnetization, double staggeredMagnetization)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(yukawa_N==0.0){ return 0.0; }
	if(!(L0==L1 && L0==L2))
	{
		std::cerr <<"Error, no L^3xL_t lattice given in constrainedEffectivePotential::computeFermionicContribution_onlyFunction_LcubeTimesLt" <<std::endl;
		exit(EXIT_FAILURE);
	}
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
// 	double one_ov_twoRho=0.5/rho;
	double ySquared(yukawa_N*yukawa_N);
	double ySquaredmSquared(ySquared*mSq);
	double ySquaredsSquared(ySquared*sSq);
	//dummies needed for the summation
	double dummyForAddition(0.0);
	int Lhalf=L0/2;
	int LT_half=L3/2;
	for(int l3=1; l3<LT_half; ++l3)
	{
		for(int l0=1; l0<Lhalf; ++l0)
		{
			for(int l1=l0+1; l1<Lhalf; ++l1)
			{
				for(int l2=l1+1; l2<Lhalf; ++l2)
				{
					//p,q,r 48
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,48.0);
				}
				{
					int l2=l1;
					//p,q,q, 24
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,24.0);
					l2=l0;
					//p,p,q, 24
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,24.0);
					l2=Lhalf;
					//p,q,L/2 24
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,24.0);
					l2=0;
					//0,p,q, 24
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,24.0);
				}
			}
			{
				int l1=l0;
				{
					int l2=l1;
					//p,p,p, 8
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,8.0);
					l2=Lhalf;
					//p,p,L/2 12
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,12.0);
					l2=0;
					//0,p,p, 12
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,12.0);
				}
				l1=Lhalf;
				{
					int l2=Lhalf;
					//p,L/2,L/2, 6
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,6.0);
					l2=0;
					//0,p,L/2, 12
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,12.0);
				}
				l1=0;
				{
					int l2=0;
					//0,0,p  6
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,6.0);
				}
			}
		}
		{
			int l0=Lhalf;
			{
				int l1=Lhalf;
				{
					int l2=Lhalf;
					//L/2,L/2,L/2, 1
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,1.0);
					l2=0;
					//0,L/2,L/2, 3
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,3.0);
				}
				l1=0;
				{
					int l2=0;
					//0,0,L/2, 3
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,3.0);
				}
			}
			l0=0;
			{
				int l1=0;
				{
					int l2=0;
					//0,0,0, 1
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,1.0);
				}
			}
		}
	}
	//now multiplzy with two, due to the symmetry in L3
	dummyForAddition*=2.0;
	//now do the same thing for l3=0 and l3=Lt_half
	double BC_factor=(antiperiodicBC_L3 ? 2.0 :1.0);
	for(int l3=0; l3<=LT_half; l3+=(antiperiodicBC_L3 ? LT_half*2 :LT_half))
	{
		for(int l0=1; l0<Lhalf; ++l0)
		{
			for(int l1=l0+1; l1<Lhalf; ++l1)
			{
				for(int l2=l1+1; l2<Lhalf; ++l2)
				{
					//p,q,r 48
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,48.0)*BC_factor;
				}
				{
					int l2=l1;
					//p,q,q, 24
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,24.0)*BC_factor;
					l2=l0;
					//p,p,q, 24
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,24.0)*BC_factor;
					l2=Lhalf;
					//p,q,L/2 24
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,24.0)*BC_factor;
					l2=0;
					//0,p,q, 24
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,24.0)*BC_factor;
				}
			}
			{
				int l1=l0;
				{
					int l2=l1;
					//p,p,p, 8
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,8.0)*BC_factor;
					l2=Lhalf;
					//p,p,L/2 12
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,12.0)*BC_factor;
					l2=0;
					//0,p,p, 12
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,12.0)*BC_factor;
				}
				l1=Lhalf;
				{
					int l2=Lhalf;
					//p,L/2,L/2, 6
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,6.0)*BC_factor;
					l2=0;
					//0,p,L/2, 12
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,12.0)*BC_factor;
				}
				l1=0;
				{
					int l2=0;
					//0,0,p  6
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,6.0)*BC_factor;
				}
			}
		}
		{
			int l0=Lhalf;
			{
				int l1=Lhalf;
				{
					int l2=Lhalf;
					//L/2,L/2,L/2, 1
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,1.0)*BC_factor;
					l2=0;
					//0,L/2,L/2, 3
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,3.0)*BC_factor;
				}
				l1=0;
				{
					int l2=0;
					//0,0,L/2, 3
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,3.0)*BC_factor;
				}
			}
			l0=0;
			{
				int l1=0;
				{
					int l2=0;
					//0,0,0, 1
					dummyForAddition+=fermionicContributionInline_onlyFunction(l0,l1,l2,l3,ySquaredmSquared,ySquaredsSquared,1.0)*BC_factor;
				}
			}
		}
	}
	
	dummyForAddition*=-2.0;
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
	return dummyForAddition;
}
*/


/*
double constrainedEffectivePotential::fermionicContributionInline_onlyFunction(int l0, int l1, int l2, int l3, double ySq_mSq, double ySq_sSq, double factor)
{
	std::complex< double > nuOfP(0.,0.), nuOfVarP(0.,0.);
	double abs_nuOfP, abs_nuOfVarP;
	double abs_gammaOfP, abs_gammaOfVarP;
	double a,b;
	
	computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
	abs_nuOfP=std::abs(nuOfP); 
	abs_nuOfVarP=std::abs(nuOfVarP);
	abs_gammaOfP=std::abs( 1.0 - one_ov_twoRho*nuOfP );
	abs_gammaOfVarP=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
	a = abs_nuOfP*abs_nuOfVarP + (ySq_mSq - ySq_sSq)*abs_gammaOfP*abs_gammaOfVarP;
	b = abs_gammaOfP*abs_nuOfVarP - abs_nuOfP*abs_gammaOfVarP;
	return factor*log(a*a + ySq_mSq*b*b);
}
*/

