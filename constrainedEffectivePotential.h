#ifndef CONSTRAINEDEFFECTIVEPOTENTIAL_H
#define CONSTRAINEDEFFECTIVEPOTENTIAL_H


// implements the computation of ter CEP as in Philipp's thesis chapter 4.
// there will be the possibility to implement antiperiodic boundary conditions in time 
// for the L3 direction
// as in the thesis, only mass degenerate quarks

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "gsl/gsl_multimin.h"


class constrainedEffectivePotential
{
	private:
	//physical parameters
	int L0, L1, L2, L3; //only L3 can be chosen antiperiodic
	bool antiperiodicBC_L3;
	
	double *sinSquaredOfPmu_L0, *sinSquaredOfPmu_L1, *sinSquaredOfPmu_L2, *sinSquaredOfPmu_L3;
	double *sinSquaredOfPmuHalf_L0, *sinSquaredOfPmuHalf_L1, *sinSquaredOfPmuHalf_L2, *sinSquaredOfPmuHalf_L3;
	double *cosSquaredOfPmuHalf_L0, *cosSquaredOfPmuHalf_L1, *cosSquaredOfPmuHalf_L2, *cosSquaredOfPmuHalf_L3;
	
	double kappa_N; //==kappa
	double lambda_N; //==lambda * N_f
	double yukawa_N; //== Y_N= Y*\sqrt(N_f)
	
	int N_f;
	double rho;
	double r;
	
	//parameters for the minimization
	double tolerance;
	int minimizationAlgorithm;
	int maxNumerOfIterations;
	
	//gsl stuff for minimization
	// a minimizer using the function and it's gradient
	gsl_multimin_fdfminimizer *minimizer;
	gsl_multimin_function_fdf functionHandler;
	
	public:
	//constructor
	constrainedEffectivePotential(int L0, int L1, int L2, int L3, bool antiL3=false);
	//destructor
	~constrainedEffectivePotential();
	
	//setting values
	void set_kappa_N(double new_k);
	void set_lambda_N(double new_l);
	void set_yukawa_N(double new_y);
	
	void set_N_f(int N);
	void set_rho(double new_rho);
	void set_r(double new_r);
	
	void set_tolerance(double new_tol);
	void set_minimizationAlgorithm(int new_alg);
	void set_maxNumerOfIterations(int new_NOI);
	
	
	//private functions
// 	private: NOTE make private in the end
	void fillLatticeMomenta();

	//computation of eigenvalues
	std::complex< double > computeAnalyticalEigenvalue(double p0, double p1, double p2, double p3);
	//computes eigenvalues from the index and uses the arrays to not compute the sins
	std::complex< double > computeAnalyticalEigenvalue_fromIndex(int l0, int l1, int l2, int l3);
	//computes at the same time nu(p) and nu(p+(pi,pi,pi,pi)) from the index
	void computeAnalyticalEigenvalue_fromIndex_pAndVarP(int l0, int l1, int l2, int l3, std::complex< double > &nuOfP, std::complex< double > &nuOfVarP);
	
	
	double computeConstrainedEffectivePotential_onlyFunction(double magnetization, double staggeredMagnetization);
	//versions of computing the function only
	//computation of the fermionic contribution
	//just computes the sum of logs as in philipp's thesis 4.35
	double computeFermionicContribution_onlyFunction_qad(double magnetization, double staggeredMagnetization);
	
	void computeConstrainedEffectivePotential_onlyGradient( double magnetization, double staggeredMagnetization, double &dU_ov_dm, double &dU_ov_ds);
	//versions of computing the gradient only
	void computeFermionicContribution_onlyGradient_qad( double magnetization, double staggeredMagnetization, double &dUf_ov_dm, double &dUf_ov_ds);
	
// 	void computeConstrainedEffectivePotential_FunctionAndGradient( double magnetization, double staggeredMagnetization, double &dU_ov_dm, double &du_ov_ds);
	//versions of computing the gradient only
// 	void computeConstrainedEffectivePotential_FunctionAndGradient_qad( double magnetization, double staggeredMagnetization, double &dU_ov_dm, double &du_ov_ds);
	
	

// 	void computeConstrainedEffectivePotential_gradient_qad( double magnetization, double staggeredMagnetization, double &result_magnetization, double &result_staggered );

// 	void computeConstrainedEffectivePotential_and_gradient_qad( double magnetization, double staggeredMagnetization, double &resultCEP, double &result_magnetization, double &result_staggered );
	
	
	
	
};




#endif




