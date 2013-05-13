#include "constrainedEffectivePotential.h"


double constrainedEffectivePotential::computeConstrainedEffectivePotential_onlyFunction(const double magnetization, const double staggeredMagnetization)
{
	//eq 4.35 philipp's thesis
	//U= -8*kappa*(m^2-s^2) + m^2 + s^2 + lambda*( m^4 + s^4 + 6*m^2s^2 - 2*(m^2+s^2)) + U_f(m,s)
	double result(0.0);
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
	
	//NOTE optimize here
	double fermionicContribution_sumOfLogs=computeFermionicContribution_onlyFunction_qad(magnetization,staggeredMagnetization);
	result+=fermionicContribution_sumOfLogs;
	result+= mSq + sSq - 8.0 * kappa_N * (mSq - sSq);
	result+= lambda_N * ( mSq*mSq + sSq*sSq + 6.0*mSq*sSq -2.0*(mSq + sSq) );
	
	return result;
}



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



double constrainedEffectivePotential::fermionicContributionInline_onlyFunction_FromStoredEigenvalues(int index, double ySq_mSq, double ySq_sSq)
{
	
	double a,b;
	a = absNuP[index]*absNuVarP[index] + (ySq_mSq - ySq_sSq)*absGammaP[index]*absGammaVarP[index];
	b = absGammaP[index]*absNuVarP[index] - absNuP[index]*absGammaVarP[index];
	return factorOfMomentum[index]*log(a*a + ySq_mSq*b*b);
}


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


double constrainedEffectivePotential::computeConstrainedEffectivePotential_onlyFunction_gsl(const gsl_vector *mags, void *params)
{
	double magnetization=gsl_vector_get(mags,0);
	double staggeredMagnetization=gsl_vector_get(mags,1);
	
	return computeConstrainedEffectivePotential_onlyFunction(magnetization, staggeredMagnetization);
}

double wrapper_computeConstrainedEffectivePotential_onlyFunction_gsl(const gsl_vector *mags, void *params)
{
	constrainedEffectivePotential *CEP = (constrainedEffectivePotential *)params;
	return CEP->computeConstrainedEffectivePotential_onlyFunction_gsl(mags, NULL);
}


