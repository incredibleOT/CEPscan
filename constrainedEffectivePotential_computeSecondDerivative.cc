#include "constrainedEffectivePotential.h"


void constrainedEffectivePotential::computeConstrainedEffectivePotential_secondDerivatives(const double magnetization, const double staggeredMagnetization, double &d2U_ov_dmdm, double &d2U_ov_dsds, double &d2U_ov_dmds)
{
	//U= -8*kappa*(m^2-s^2) + m^2 + s^2 + lambda*( m^4 + s^4 + 6*m^2s^2 - 2*(m^2+s^2)) + lambda_6(m^6 + s^6 + 15*m^4*s^2 + 15*m^2*s^4) + U_f(m,s)
	//dU/dm = dU_f/d_m - 16 kappa * m + 2*m + lambda*(4*m^3 + 12 *m*s^2 - 4*m) + lambda_6*(6*m^5 + 60*m^3*s^2 + 30*m*s^4)
	//dU/ds = dU_f/d_m + 16 kappa * s + 2*s + lambda*(4*s^3 + 12 *m^2*s - 4*s) + lambda_6*(6*s^5 + 30*m^4*s + 60*m^2*s^3)
	//
	//d^2U/dm^2 = d^2U_f/d_m^2 - 16 kappa + 2 + lambda*(12*m^2 + 12*s^2 - 4) + lambda_6*(30*m^4 + 180*m^2*s^2 + 30*s^4)
	//d^2U/ds^2 = d^2U_f/d_m^2 + 16 kappa + 2 + lambda*(12*s^2 + 12*m^2 - 4) + lambda_6*(30*s^4 + 180*m^2*s^2 + 30*m^4)
	//d^2U/dmds = d^2U_f/dmds  + lambda*(24*m*s) + lambda_6*(120*m^3*s + 120*m*s^3)
	//
	//NOTE added a term coming from the boson loop:
	//U+=16*(m^2+s^2)*lambda_N*N_f^{-1}*P_B
	//P_B= 1/V * sum_{p} 1/(2-4*lambda_N-4*kappa*sum{cos(P_mu)}) excluding zero and staggered mode
	//With that:
	//dU/dm +=32*m*lambda_N*N_f^{-1}*P_B
	//dU/ds +=32*s*lambda_N*N_f^{-1}*P_B
	//d^2U/dm^2 += 32*lambda_N*N_f^{-1}*P_B
	//d^2U/ds^2 += 32*lambda_N*N_f^{-1}*P_B
	//d^2U/dmds+=0;
	d2U_ov_dmdm=0.0; d2U_ov_dsds=0.0; d2U_ov_dmds=0.0;
	if(useBosonicLoop && !bosonicLoopSet){ bosonicLoop=computeBosonicPropagatorSum_fromStoredSumOfCos();}
	
	compute_fermionicContribution_secondDerivatives_FromStoredEigenvalues(magnetization, staggeredMagnetization, d2U_ov_dmdm, d2U_ov_dsds, d2U_ov_dmds);
	double mSq(magnetization*magnetization), sSq(staggeredMagnetization*staggeredMagnetization);
	d2U_ov_dmdm+= -16.0*kappa_N + 2.0;
	d2U_ov_dmdm+= lambda_N*( 12.0*(mSq + sSq) - 4.0 );
	d2U_ov_dmdm+= lambda_6_N*( 30.0*mSq*mSq + 180.0*mSq*sSq + 30.0*sSq*sSq);
	if(useBosonicLoop){ d2U_ov_dmdm += 32.0*lambda_N*bosonicLoop/static_cast< double >(N_f); }
	
	d2U_ov_dsds+= +16.0*kappa_N + 2.0;
	d2U_ov_dsds+= lambda_N*( 12.0*(mSq + sSq) - 4.0 );
	d2U_ov_dsds+= lambda_6_N*( 30.0*sSq*sSq + 180.0*mSq*sSq + 30.0*mSq*mSq);
	if(useBosonicLoop){ d2U_ov_dsds += 32.0*lambda_N*bosonicLoop/static_cast< double >(N_f); }
	
	d2U_ov_dmds+= 24.0*lambda_N*magnetization*staggeredMagnetization;
	d2U_ov_dmds+= lambda_6_N*(120.0*mSq*magnetization*staggeredMagnetization + 120.0*magnetization*sSq*staggeredMagnetization);
}


void constrainedEffectivePotential::fermionicContributionInline_secondDerivatives_FromStoredEigenvalues(int index, const double ySq_mSq, const double ySq_sSq, double &d2Uf_ov_dmdm, double &d2Uf_ov_dsds, double &d2Uf_ov_dmds)
{
	double a(absNuP[index]*absNuVarP[index] + (ySq_mSq - ySq_sSq)*absGammaP[index]*absGammaVarP[index]);
	double b(absGammaP[index]*absNuVarP[index] - absNuP[index]*absGammaVarP[index]);
	double one_ov_A(1.0/(a*a + ySq_mSq*b*b));
	//std::cout <<"a=" <<a <<"   b=" <<b <<"   one_ov_A =" <<one_ov_A <<std::endl;
	double gg(absGammaP[index]*absGammaVarP[index]);
	double two_gg_a_pl_bSq(2.0*gg*a + b*b); 
	
	d2Uf_ov_dmdm=one_ov_A*(two_gg_a_pl_bSq + 2.0*ySq_mSq*( 2.0*gg*gg - two_gg_a_pl_bSq*two_gg_a_pl_bSq*one_ov_A) )*factorOfMomentum[index];
	d2Uf_ov_dsds=gg*one_ov_A*(a + 2.0*ySq_sSq*gg*( 2.0*a*a*one_ov_A - 1.0))*factorOfMomentum[index];
	d2Uf_ov_dmds=gg*one_ov_A*( gg - a*two_gg_a_pl_bSq*one_ov_A )*factorOfMomentum[index];
// 	asdas
	//NOTE change here
}

void constrainedEffectivePotential::compute_fermionicContribution_secondDerivatives_FromStoredEigenvalues(const double magnetization, const double staggeredMagnetization, double &d2Uf_ov_dmdm, double &d2Uf_ov_dsds, double &d2Uf_ov_dmds)
{
	//U_f = -2/V*sum log[ a^2 + y^2*m^2*b^2 ] = -2/v*sum log[A]
	//a = a(m,s) = |nu(p)|*|nu(varP)| + y^2*(m^2-s^2)*|gamma(p)|*|gamma(varP)|
	//b = |gamma(p)|*|nu(varP)| - |nu(p)|*|gamma(varP)|
	//A = a^2 + y^2*m^2*b^2
	//dU_f/dm = -2/V*sum[ (4*y^2*m * |gamma(p)|*|gamma(varP)|*a + 2*m*y^2*b^2)/(a^2 + y^2*m^2*b^2 ) ] 
	//        = -4*y^2*m/V * sum[ (2* |gamma(p)|*|gamma(varP)|*a + b^2)/(a^2 + y^2*m^2*b^2 ) ]
	//        = -4*y^2*m/V * sum[ (2* |gamma(p)|*|gamma(varP)|*a + b^2)/A ]
	//dU_f/ds = -2/V*sum[ (-4*y^2*s * |gamma(p)|*|gamma(varP)|*a )/(a^2 + y^2*m^2*b^2 ) ]
	//        = +8*y^2*s/V * sum[ ( |gamma(p)|*|gamma(varP)|*a )/(a^2 + y^2*m^2*b^2 ) ]
	//        = 8*y^2*m/V * sum[ ( |gamma(p)|*|gamma(varP)*a|)/A ]
	//
	// (gg) = |gamma(p)|*|gamma(varP)|
	//d^2U_f/dmdm = -4*y^2/V * sum[ ( 2(gg)a+b^2 + 2*y^2*m^2*( 2(gg)^2 - (2(gg)a + b^2)^2/A ) )/A ]
	//
	//d^2U_f/dsds = 8*y^2/V * sum[ (gg)/A * ( a + 2*y^2*s^2*(gg)*(-1 + 2*a^2/A))]
	//
	//d^2U/dmds = 16y^4*m*s/V * sum[ (gg)/A * ( gg - a*(2(gg)a + b^2)/A)]
	if(yukawa_N < 0.0)
	{
		std::cerr <<"Error, no yukawa coupling set in constrainedEffectivePotential::compute_fermionicContribution_secondDerivatives_FromStoredEigenvalues(const double magnetization, const double staggeredMagnetization, double &d2Uf_ov_dmdm, double &d2Uf_ov_dsds, double &d2Uf_ov_dmds)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	d2Uf_ov_dmdm=0.0; d2Uf_ov_dsds=0.0; d2Uf_ov_dmds=0.0;
	if(yukawa_N==0.0){ return; }
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
	double ySquared(yukawa_N*yukawa_N);
	double ySquaredmSquared(ySquared*mSq);
	double ySquaredsSquared(ySquared*sSq);
	double dummyForAddition_dmdm(0.0), dummyForAddition_dsds(0.0), dummyForAddition_dmds(0.0);
	for(int index=0; index<numberOfDistingtMomenta; ++index)
	{
		fermionicContributionInline_secondDerivatives_FromStoredEigenvalues(index, ySquaredmSquared, ySquaredsSquared, dummyForAddition_dmdm, dummyForAddition_dsds, dummyForAddition_dmds);
		d2Uf_ov_dmdm+=dummyForAddition_dmdm;
		d2Uf_ov_dsds+=dummyForAddition_dsds;
		d2Uf_ov_dmds+=dummyForAddition_dmds;
	}
	d2Uf_ov_dmdm*=-4.0*ySquared;
	d2Uf_ov_dmdm/=static_cast< double >(L0); d2Uf_ov_dmdm/=static_cast< double >(L1);
	d2Uf_ov_dmdm/=static_cast< double >(L2); d2Uf_ov_dmdm/=static_cast< double >(L3);
	
	d2Uf_ov_dsds*=8.0*ySquared;
	d2Uf_ov_dsds/=static_cast< double >(L0); d2Uf_ov_dsds/=static_cast< double >(L1);
	d2Uf_ov_dsds/=static_cast< double >(L2); d2Uf_ov_dsds/=static_cast< double >(L3);
	
	d2Uf_ov_dmds*=16.0*ySquared*ySquared*magnetization*staggeredMagnetization;
	d2Uf_ov_dmds/=static_cast< double >(L0); d2Uf_ov_dmds/=static_cast< double >(L1);
	d2Uf_ov_dmds/=static_cast< double >(L2); d2Uf_ov_dmds/=static_cast< double >(L3);
}





