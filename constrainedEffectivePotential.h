#ifndef CONSTRAINEDEFFECTIVEPOTENTIAL_H
#define CONSTRAINEDEFFECTIVEPOTENTIAL_H


// implements the computation of ter CEP as in Philipp's thesis chapter 4.
// there will be the possibility to implement antiperiodic boundary conditions in time 
// for the L3 direction
// as in the thesis, only mass degenerate quarks


#include "gsl_multimin.h"


class constrainedEffectivePotential
{
	private:
	//physical parameters
	int L0, L1, L2, L3; //only L3 can be chosen antiperiodic
	bool antiperiodicBC_L3;
	
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
	
	
	//private functions
	private:
	double computeConstrainedEffectivePotential_qad( double magnetization, double staggeredMagnetization );
	void computeConstrainedEffectivePotential_gradient_qad( double magnetization, double staggeredMagnetization, double &result_magnetization, double &result_staggered );
	void computeConstrainedEffectivePotential_and_gradient_qad( double magnetization, double staggeredMagnetization, double &resultCEP, double &result_magnetization, double &result_staggered );
	
	
	
};




#endif




