#include "constrainedEffectivePotential.h"




double constrainedEffectivePotential::computeConstrainedEffectivePotential_onlyFunction(double magnetization, double staggeredMagnetization)
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

//computation of the fermionic contribution
//just computes the sum of logs as in philipp's thesis 4.35
double constrainedEffectivePotential::computeFermionicContribution_onlyFunction_qad(double magnetization, double staggeredMagnetization)
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
	double one_ov_twoRho=0.5/rho;
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
					//std::cout <<"indexEVcomputation:   nu(p)=" <<nuOfP <<"   abs(nu(p))=" << abs_nuOfP
					//<<"   nu(var(p))="<<nuOfVarP <<"   abs(nu(varP))=" << abs_nuOfVarP <<std::endl;
					a = abs_nuOfP*abs_nuOfVarP + ySquared*(mSq - sSq)*abs_gammaOfP*abs_gammaOfVarP;
					b = abs_gammaOfP*abs_nuOfVarP + abs_nuOfP*abs_gammaOfVarP;
					dummyForAddition+=log(a*a + ySquared*mSq*b*b);
					//std::cout <<"dummyForAddition=" <<dummyForAddition <<std::endl;
				}
			}
		}
	}
	dummyForAddition*=-2.0;
// 	dummyForMultiplication=-2.0*log(dummyForMultiplication);
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
// 	dummyForMultiplication/=static_cast< double >(L0); dummyForMultiplication/=static_cast< double >(L1); dummyForMultiplication/=static_cast< double >(L2); dummyForMultiplication/=static_cast< double >(L3); 
	//std::cout <<"U_f from addition of logs = " <<dummyForAddition <<std::endl;
	return dummyForAddition;
}
