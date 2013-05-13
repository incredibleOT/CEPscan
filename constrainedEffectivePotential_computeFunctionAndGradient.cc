#include "constrainedEffectivePotential.h"


void constrainedEffectivePotential::computeConstrainedEffectivePotential_FunctionAndGradient(double const magnetization, const double staggeredMagnetization, double &U, double &dU_ov_dm, double &dU_ov_ds)
{
	//U= -8*kappa*(m^2-s^2) + m^2 + s^2 + lambda*( m^4 + s^4 + 6*m^2s^2 - 2*(m^2+s^2)) + U_f(m,s)
	//dU/dm = dU_f/d_m - 16 kappa * m + 2*m + lambda*(4*m^3 + 12 *m*s^2 - 4*m)
	//dU/ds = dU_f/d_m + 16 kappa * s + 2*s + lambda*(4*s^3 + 12 *m^2*s - 4*s)
	U=0.0; dU_ov_dm=0.0; dU_ov_ds=0.0;
	computeFermionicContribution_FunctionAndGradient_qad(magnetization, staggeredMagnetization, U, dU_ov_dm, dU_ov_ds);
	double mSq=magnetization*magnetization;
	double sSq=staggeredMagnetization*staggeredMagnetization;
	
	//NOTE optimize here
	U+= mSq + sSq - 8.0 * kappa_N * (mSq - sSq);
	U+= lambda_N * ( mSq*mSq + sSq*sSq + 6.0*mSq*sSq -2.0*(mSq + sSq) );
	
	dU_ov_dm += -16.0*kappa_N*magnetization + 2.0*magnetization;
	dU_ov_dm += lambda_N*( 4.0*magnetization*magnetization*magnetization + 12.0*magnetization*staggeredMagnetization*staggeredMagnetization - 4.0*magnetization);
	
	dU_ov_ds += +16.0*kappa_N*staggeredMagnetization + 2.0*staggeredMagnetization;
	dU_ov_ds += lambda_N*( 4.0*staggeredMagnetization*staggeredMagnetization*staggeredMagnetization + 12.0*magnetization*magnetization*staggeredMagnetization - 4.0*staggeredMagnetization);
	
	
}


void constrainedEffectivePotential::computeFermionicContribution_FunctionAndGradient_qad(const double magnetization, const double staggeredMagnetization, double &Uf, double &dUf_ov_dm, double &dUf_ov_ds)
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
		std::cerr <<"Error, no yukawa coupling set in constrainedEffectivePotential::computeFermionicContribution_FunctionAndGradient_qad(double magnetization, double staggeredMagnetization, double &Uf, double &dUf_ov_dm, double &dUf_ov_ds)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	Uf=0.0; dUf_ov_dm=0.0; dUf_ov_ds=0.0; 
	if(yukawa_N==0.0){ return; }
	
// 	double one_ov_twoRho=0.5/rho;
	
	double ySquared(yukawa_N*yukawa_N);
	double ySq_times_mSq_mi_sSq( ySquared*(magnetization*magnetization - staggeredMagnetization*staggeredMagnetization) );
	double ySq_times_mSq( ySquared*magnetization*magnetization );
	
	std::complex< double > nuOfP(0.,0.), nuOfVarP(0.,0.);
	double abs_nuOfP, abs_nuOfVarP;
	double abs_gammaOfP, abs_gammaOfVarP;
	double a,b,A,one_ov_A;
	
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
					A=(a*a + ySq_times_mSq*b*b);
					if(A==0.0)
					{
						std::cerr <<"Error, argument of log is zero for:" <<std::endl;
						std::cerr <<"l0=" <<l0 <<"  l1=" <<l1 <<"  l2=" <<l2 <<"  l3=" <<l3 <<std::endl;
						std::cerr <<"nu(p)=" <<nuOfP <<"  nu(varP)=" <<nuOfVarP <<"  gamma(p)=" <<( 1.0 - one_ov_twoRho*nuOfP ) <<"  gamma(varP)=" <<( 1.0 - one_ov_twoRho*nuOfVarP ) <<std::endl;
					}
					one_ov_A=1.0/A;
					Uf+=log(A);
					dUf_ov_dm+=(2.0*abs_gammaOfP*abs_gammaOfVarP*a + b*b)*one_ov_A;
					dUf_ov_ds+=(abs_gammaOfP*abs_gammaOfVarP*a)*one_ov_A;
				}
			}
		}
	}
	Uf*=-2.0;
	dUf_ov_dm*=-4.0*ySquared*magnetization;
	dUf_ov_ds*=8.0*ySquared*staggeredMagnetization;
	
	Uf/=static_cast< double >(L0); Uf/=static_cast< double >(L1); Uf/=static_cast< double >(L2); Uf/=static_cast< double >(L3);
	dUf_ov_dm/=static_cast< double >(L0); dUf_ov_dm/=static_cast< double >(L1); dUf_ov_dm/=static_cast< double >(L2); dUf_ov_dm/=static_cast< double >(L3);
	dUf_ov_ds/=static_cast< double >(L0); dUf_ov_ds/=static_cast< double >(L1); dUf_ov_ds/=static_cast< double >(L2); dUf_ov_ds/=static_cast< double >(L3);
	
}

void constrainedEffectivePotential::computeConstrainedEffectivePotential_FunctionAndGradient_gsl(const gsl_vector *mags, void *params, double *U, gsl_vector *gradient_of_U)
{
	double magnetization=gsl_vector_get(mags,0);
	double staggeredMagnetization=gsl_vector_get(mags,1);
	
	double dU_dm(0.0), dU_ds(0.0);
	
	computeConstrainedEffectivePotential_FunctionAndGradient(magnetization, staggeredMagnetization, *U, dU_dm, dU_ds);
	gsl_vector_set(gradient_of_U, 0, dU_dm);
	gsl_vector_set(gradient_of_U, 1, dU_ds);
}

void wrapper_computeConstrainedEffectivePotential_FunctionAndGradient_gsl(const gsl_vector *mags, void *params, double *U, gsl_vector *gradient_of_U)
{
	constrainedEffectivePotential *CEP = (constrainedEffectivePotential *)params;
	CEP->computeConstrainedEffectivePotential_FunctionAndGradient_gsl(mags, NULL, U, gradient_of_U);
}


