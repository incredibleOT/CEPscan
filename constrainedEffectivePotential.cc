#include "constrainedEffectivePotential.h"

#include "constrainedEffectivePotential_computeFunctionOnly.cc"
#include "constrainedEffectivePotential_computeGradientOnly.cc"
#include "constrainedEffectivePotential_FunctionAndGradient.cc"

constrainedEffectivePotential::constrainedEffectivePotential(int l0, int l1, int l2, int l3, bool antiL3):
L0(l0), L1(l1), L2(l2), L3(l3), antiperiodicBC_L3(antiL3),
sinSquaredOfPmu_L0(NULL),sinSquaredOfPmu_L1(NULL),sinSquaredOfPmu_L2(NULL),sinSquaredOfPmu_L3(NULL),
sinSquaredOfPmuHalf_L0(NULL),sinSquaredOfPmuHalf_L1(NULL),sinSquaredOfPmuHalf_L2(NULL),sinSquaredOfPmuHalf_L3(NULL),
cosSquaredOfPmuHalf_L0(NULL),cosSquaredOfPmuHalf_L1(NULL),cosSquaredOfPmuHalf_L2(NULL),cosSquaredOfPmuHalf_L3(NULL),
kappa_N(-1.0), lambda_N(-1.0), yukawa_N(-1.0),
N_f(1), rho(1.0), r(0.5),
tolerance(-1.0), minimizationAlgorithm(-1), maxNumerOfIterations(-1),
minimizer(NULL)
{
	if( L0%2!=0 || L1%2!=0 || L2%2!=0 || L3%2!=0 )
	{
		std::cerr <<"Error, only even Lattice extends possible!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	minimizer = new gsl_multimin_fdfminimizer;
	fillLatticeMomenta();
}


constrainedEffectivePotential::~constrainedEffectivePotential()
{
	delete minimizer;
	delete [] sinSquaredOfPmu_L0;
	delete [] sinSquaredOfPmu_L1;
	delete [] sinSquaredOfPmu_L2;
	delete [] sinSquaredOfPmu_L3;
	delete [] sinSquaredOfPmuHalf_L0;
	delete [] sinSquaredOfPmuHalf_L1;
	delete [] sinSquaredOfPmuHalf_L2;
	delete [] sinSquaredOfPmuHalf_L3;
	delete [] cosSquaredOfPmuHalf_L0;
	delete [] cosSquaredOfPmuHalf_L1;
	delete [] cosSquaredOfPmuHalf_L2;
	delete [] cosSquaredOfPmuHalf_L3;
}

void constrainedEffectivePotential::fillLatticeMomenta()
{
	double one_ov_L0(1.0/L0);
	double one_ov_L1(1.0/L1);
	double one_ov_L2(1.0/L2);
	double one_ov_L3(1.0/L3);
	const double PI(atan(1) * 4.0);
	double toAddForL3 = antiperiodicBC_L3 ? PI/L3 : 0.0;
	sinSquaredOfPmu_L0 = new double [L0]; sinSquaredOfPmuHalf_L0 = new double [L0]; cosSquaredOfPmuHalf_L0 = new double [L0];
	sinSquaredOfPmu_L1 = new double [L1]; sinSquaredOfPmuHalf_L1 = new double [L1]; cosSquaredOfPmuHalf_L1 = new double [L1];
	sinSquaredOfPmu_L2 = new double [L2]; sinSquaredOfPmuHalf_L2 = new double [L2]; cosSquaredOfPmuHalf_L2 = new double [L2];
	sinSquaredOfPmu_L3 = new double [L3]; sinSquaredOfPmuHalf_L3 = new double [L3]; cosSquaredOfPmuHalf_L3 = new double [L3];
	double p,sinP,cosP;
	for(int i=0; i<L0; ++i)
	{
		p=2.0 * PI * i * one_ov_L0; 
		sinP=sin(p); sinSquaredOfPmu_L0[i]=sinP*sinP;
		sinP=sin(0.5*p); sinSquaredOfPmuHalf_L0[i]=sinP*sinP;
		cosP=cos(0.5*p); cosSquaredOfPmuHalf_L0[i]=cosP*cosP;
		std::cout <<"l0 = " << i <<"   p_mu = " <<p <<"   sinSquaredOfPmu = " <<sinSquaredOfPmu_L0[i] <<"   sinSquaredOfPmuHalf = " <<sinSquaredOfPmuHalf_L0[i] <<std::endl;
	}
	for(int i=0; i<L1; ++i)
	{
		p=2.0 * PI * i * one_ov_L1; 
		sinP=sin(p); sinSquaredOfPmu_L1[i]=sinP*sinP;
		sinP=sin(0.5*p); sinSquaredOfPmuHalf_L1[i]=sinP*sinP;
		cosP=cos(0.5*p); cosSquaredOfPmuHalf_L1[i]=cosP*cosP;
		std::cout <<"l1 = " << i <<"   p_mu = " <<p <<"   sinSquaredOfPmu = " <<sinSquaredOfPmu_L1[i] <<"   sinSquaredOfPmuHalf = " <<sinSquaredOfPmuHalf_L1[i] <<std::endl;
	}
	for(int i=0; i<L2; ++i)
	{
		p=2.0 * PI * i * one_ov_L2; 
		sinP=sin(p); sinSquaredOfPmu_L2[i]=sinP*sinP;
		sinP=sin(0.5*p); sinSquaredOfPmuHalf_L2[i]=sinP*sinP;
		cosP=cos(0.5*p); cosSquaredOfPmuHalf_L2[i]=cosP*cosP;
		std::cout <<"l2 = " << i <<"   p_mu = " <<p <<"   sinSquaredOfPmu = " <<sinSquaredOfPmu_L2[i] <<"   sinSquaredOfPmuHalf = " <<sinSquaredOfPmuHalf_L2[i] <<std::endl;
	}
	for(int i=0; i<L3; ++i)
	{
		p=2.0 * PI * i * one_ov_L3 + toAddForL3; 
		sinP=sin(p); sinSquaredOfPmu_L3[i]=sinP*sinP;
		sinP=sin(0.5*p); sinSquaredOfPmuHalf_L3[i]=sinP*sinP;
		cosP=cos(0.5*p); cosSquaredOfPmuHalf_L3[i]=cosP*cosP;
		std::cout <<"l3 = " << i <<"   p_mu = " <<p <<"   sinSquaredOfPmu = " <<sinSquaredOfPmu_L3[i] <<"   sinSquaredOfPmuHalf = " <<sinSquaredOfPmuHalf_L3[i] <<std::endl;
	}
}

void constrainedEffectivePotential::set_kappa_N(double new_k){ kappa_N=new_k; }
void constrainedEffectivePotential::set_lambda_N(double new_l){ lambda_N=new_l; }
void constrainedEffectivePotential::set_yukawa_N(double new_y){ yukawa_N=new_y; }

void constrainedEffectivePotential::set_N_f(int new_N){ N_f=new_N; }
void constrainedEffectivePotential::set_rho(double new_rho){ rho=new_rho; }
void constrainedEffectivePotential::set_r(double new_r){ r=new_r; }

void constrainedEffectivePotential::set_tolerance(double new_tol){ tolerance=new_tol; }
void constrainedEffectivePotential::set_minimizationAlgorithm(int new_alg){ minimizationAlgorithm=new_alg; }
void constrainedEffectivePotential::set_maxNumerOfIterations(int new_NOI){ maxNumerOfIterations=new_NOI; }





std::complex< double > constrainedEffectivePotential::computeAnalyticalEigenvalue(double p0, double p1, double p2, double p3)
{
	//computes \nu^{+} from philipp's thesis (eq 3.9)
	// \nu(p) = \rho/a + \rho/a * [i\sqrt{\tilde{p}^2} + a*r *\hat{p}^2 - \rho/a]/[\srqt{\tilde{p}^2 + (a*r *\hat{p}^2 - \rho/a)^2}]
	//a is set to 1
	//\hat{p}^2 = 1/a^2 \sum_{\mu} 4 \sin{a p_\mu/2}^2
	//\tilde{p}^2 = 1/a^2 \sum_{\mu} \sin{a p_\mu}^2
	//NOTE: in philipp's thesis it's r/2 NOTE May be an ambiguity in the definition of r...
	double dummy=sin(p0); double p_tilde_sq=dummy*dummy; 
	dummy=sin(p1); p_tilde_sq+=dummy*dummy;
	dummy=sin(p2); p_tilde_sq+=dummy*dummy;
	dummy=sin(p3); p_tilde_sq+=dummy*dummy;
	
	dummy=sin(0.5*p0); double p_hat_sq = dummy*dummy;
	dummy=sin(0.5*p1); p_hat_sq += dummy*dummy;
	dummy=sin(0.5*p2); p_hat_sq += dummy*dummy;
	dummy=sin(0.5*p3); p_hat_sq += dummy*dummy;
	p_hat_sq*=4.0;
	
// 	std::cout <<"old method: p_tilde_sq = " <<p_tilde_sq <<"     p_hat_sq = "<<p_hat_sq <<std::endl;
// 	std::cout <<"p0=" <<p0 <<"   p1=" <<p1 <<"   p2=" <<p2 <<"   p3=" <<p3 <<std::endl;
	
	double one_ov_denom=1.0/sqrt(p_tilde_sq + (r*p_hat_sq - rho)*(r*p_hat_sq - rho) );
	
	return std::complex< double >(rho + rho*one_ov_denom*(r*p_hat_sq - rho) , rho*one_ov_denom*sqrt(p_tilde_sq) );
}



std::complex< double > constrainedEffectivePotential::computeAnalyticalEigenvalue_fromIndex(int l0, int l1, int l2, int l3)
{
	double p_tilde_sq = sinSquaredOfPmu_L0[l0] + sinSquaredOfPmu_L1[l1] + 
	                    sinSquaredOfPmu_L2[l2] + sinSquaredOfPmu_L3[l3];
	
	
	double p_hat_sq = sinSquaredOfPmuHalf_L0[l0] + sinSquaredOfPmuHalf_L1[l1] +
	                  sinSquaredOfPmuHalf_L2[l2] + sinSquaredOfPmuHalf_L3[l3];
	p_hat_sq*=4.0;

	
	double one_ov_denom=1.0/sqrt(p_tilde_sq + (r*p_hat_sq - rho)*(r*p_hat_sq - rho) );
	return std::complex< double >(rho + rho*one_ov_denom*(r*p_hat_sq - rho) , rho*one_ov_denom*sqrt(p_tilde_sq) );
}



void constrainedEffectivePotential::computeAnalyticalEigenvalue_fromIndex_pAndVarP(int l0, int l1, int l2, int l3, std::complex< double > &nuOfP, std::complex< double > &nuOfVarP)
{
// simultaniously gives the eigenvalues for p and p+(pi,pi,pi,pi)=var(p)
// uses that sin^2(x+pi)=sin^2(x) for the tilde(p) nad
// sin^2((x+pi)/2)=cos^2(x) for hat(p)
// for varp the p_hat part is different
	double p_tilde_sq = sinSquaredOfPmu_L0[l0] + sinSquaredOfPmu_L1[l1] + 
	                    sinSquaredOfPmu_L2[l2] + sinSquaredOfPmu_L3[l3];
	
	double p_hat_sq = sinSquaredOfPmuHalf_L0[l0] + sinSquaredOfPmuHalf_L1[l1] +
	                  sinSquaredOfPmuHalf_L2[l2] + sinSquaredOfPmuHalf_L3[l3];
	p_hat_sq*=4.0;
	
	double varP_hat_sq = cosSquaredOfPmuHalf_L0[l0] + cosSquaredOfPmuHalf_L1[l1] +
	                     cosSquaredOfPmuHalf_L2[l2] + cosSquaredOfPmuHalf_L3[l3];
	varP_hat_sq*=4.0;
	
	double sqrt_p_tilde_sq(sqrt(p_tilde_sq));
	
	double one_ov_denom=1.0/sqrt(p_tilde_sq + (r*p_hat_sq - rho)*(r*p_hat_sq - rho) );
	nuOfP=std::complex< double >(rho + rho*one_ov_denom*(r*p_hat_sq - rho) , rho*one_ov_denom*sqrt_p_tilde_sq );
	
	one_ov_denom=1.0/sqrt(p_tilde_sq + (r*varP_hat_sq - rho)*(r*varP_hat_sq - rho) );
	nuOfVarP=std::complex< double >(rho + rho*one_ov_denom*(r*varP_hat_sq - rho) , rho*one_ov_denom*sqrt_p_tilde_sq );
	
}


