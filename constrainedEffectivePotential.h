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
#include "gsl/gsl_vector.h"

class constrainedEffectivePotential
{
	private:
	//physical parameters
	int L0, L1, L2, L3; //only L3 can be chosen antiperiodic
	bool antiperiodicBC_L3;
	
	//fermionic momenta
	double *sinSquaredOfPmu_L0, *sinSquaredOfPmu_L1, *sinSquaredOfPmu_L2, *sinSquaredOfPmu_L3;
	double *sinSquaredOfPmuHalf_L0, *sinSquaredOfPmuHalf_L1, *sinSquaredOfPmuHalf_L2, *sinSquaredOfPmuHalf_L3;
	double *cosSquaredOfPmuHalf_L0, *cosSquaredOfPmuHalf_L1, *cosSquaredOfPmuHalf_L2, *cosSquaredOfPmuHalf_L3;
	//for the bosonic sum
	double *cosOfPmu_L0, *cosOfPmu_L1, *cosOfPmu_L2, *cosOfPmu_L3;
	
	
	int numberOfDistingtMomenta;
	double *absNuP, *absNuVarP, *absGammaP, *absGammaVarP, *factorOfMomentum;
	int numberOfDistingtMomenta_bosonic;
	double *sumOfCosOfPmu, *factorOfMomentum_bosonic;
	
	double kappa_N; //==kappa
	double lambda_N; //==\hat lambda * N_f
	double lambda_6_N; //==\hat lambda_& * N_f^2
	double yukawa_N; //== Y_N= \hat Y*\sqrt(N_f)
	
	int N_f;
	double rho;
	double one_ov_twoRho;
	double r;
	
	//stores sum_p {1/(2-4lambda_N-4kappa*sum(cos(p_mu)))}
	double bosonicLoop;
	bool bosonicLoopSet;
	bool useBosonicLoop;
	//parameters for the minimization
	double toleranceForLineMinimization;
	double toleranceForConvergence;
	double initialStepSize;
	int maxNumerOfIterations;
	int minimizationAlgorithm; 
	bool minimizerInitialized;
	bool iteratorStoppedFlag;
	// 1 = gsl_multimin_fdfminimizer_conjugate_fr
	// 2 = gsl_multimin_fdfminimizer_conjugate_pr
	// 3 = gsl_multimin_fdfminimizer_vector_bfgs2
	// 4 = gsl_multimin_fdfminimizer_vector_bfgs
	// 5 = gsl_multimin_fdfminimizer_steepest_descent
	
	
	//gsl stuff for minimization
	// a minimizer using the function and it's gradient
	gsl_multimin_fdfminimizer *minimizer;
	gsl_multimin_function_fdf functionHandler;
	const gsl_multimin_fdfminimizer_type *algorithmForMinimization;
	
	public:
	//constructor
	constrainedEffectivePotential(const int L0, const int L1, const int L2, const int L3, const bool antiL3=false);
	//destructor
	~constrainedEffectivePotential();
	
	//setting values
	void set_kappa_N(double new_k);
	void set_lambda_N(double new_l);
	void set_lambda_6_N(double new_l_6);
	void set_yukawa_N(double new_y);
	void set_kappa_lambda_yukawa_N(double new_k, double new_l, double new_y);
	void set_kappa_lambda_lambda_6_yukawa_N(double new_k, double new_l, double new_l_6, double new_y);
	
	void set_N_f(int N);
	void set_rho(double new_rho);
	void set_r(double new_r);
	
	void set_useBosonicLoop(bool newSet);
	
	void set_toleranceForLineMinimization(double new_tol);
	void set_toleranceForConvergence(double new_tol);
	void set_initialStepSize(double new_step);
	void set_maxNumerOfIterations(int new_NOI);
	void set_minimizationAlgorithm(int new_alg);
	
	double get_kappa_N();
	double get_lambda_N();
	double get_lambda_6_N();
	double get_yukawa_N();
	
	int get_N_f();
	double get_rho();
	double get_r();
	
	double  get_toleranceForLineMinimization();
	double  get_toleranceForConvergence();
	double  get_initialStepSize();
	int  get_maxNumerOfIterations();
	int  get_minimizationAlgorithm();
	
	
	//private functions
	private: 
	void fillLatticeMomenta();
	void fillEigenvalues();
	void fillBosonicLoopValues();

	//computation of eigenvalues
	std::complex< double > computeAnalyticalEigenvalue(const double p0, const double p1, const double p2, const double p3);
	//computes eigenvalues from the index and uses the arrays to not compute the sins
	std::complex< double > computeAnalyticalEigenvalue_fromIndex(const int l0, const int l1, const int l2, const int l3);
	//computes at the same time nu(p) and nu(p+(pi,pi,pi,pi)) from the index
	void computeAnalyticalEigenvalue_fromIndex_pAndVarP(const int l0, const int l1, const int l2, const int l3, std::complex< double > &nuOfP, std::complex< double > &nuOfVarP);
	
	
	public:
	
	double computeBosonicPropagatorSum_qad();
	double computeBosonicPropagatorSum_fromStoredSumOfCos();
	
	double computeConstrainedEffectivePotential_onlyFunction(const double magnetization, const double staggeredMagnetization);
	//versions of computing the function only
	//computation of the fermionic contribution
	//just computes the sum of logs as in philipp's thesis 4.35
	double computeFermionicContribution_onlyFunction_qad(const double magnetization, const double staggeredMagnetization);
	
	//specialized for L^3xL_t
	/*
	double computeFermionicContribution_onlyFunction_LcubeTimesLt(const double magnetization, const double staggeredMagnetization);
	double computeFermionicContribution_onlyFunction_LcubeTimesLt_withInline(const double magnetization, const double staggeredMagnetization);
	double fermionicContributionInline_onlyFunction(int l0, int l1, int l2, int l3, double ySq_mSq, double ySq_sSq, double factor);
	*/ //not used 
	//uses stored values NOTE, this will be used in the end
	double computeFermionicContribution_onlyFunction_FromStoredEigenvalues(const double magnetization, const double staggeredMagnetization);
	// format gsl likes
	double computeConstrainedEffectivePotential_onlyFunction_gsl(const gsl_vector *mags);
	//for nicer code
	
	double fermionicContributionInline_onlyFunction_FromStoredEigenvalues(int index, double ySq_mSq, double ySq_sSq);
	//wrapper for function to make a pointer below the class

	
	//computes only the gradient
	void computeConstrainedEffectivePotential_onlyGradient( const double magnetization, const double staggeredMagnetization, double &dU_ov_dm, double &dU_ov_ds);
	//quick and dirty
	void computeFermionicContribution_onlyGradient_qad( const double magnetization, const double staggeredMagnetization, double &dUf_ov_dm, double &dUf_ov_ds);
	//faster
	void computeFermionicContribution_onlyGradient_FromStoredEigenvalues( const double magnetization, const double staggeredMagnetization, double &dUf_ov_dm, double &dUf_ov_ds);
	// format gsl likes
	void computeConstrainedEffectivePotential_onlyGradient_gsl(const gsl_vector *mags, gsl_vector *gradient_of_U);
	//for nicer code
	void fermionicContributionInline_onlyGradient_FromStoredEigenvalues(int index, double ySq_mSq, double ySq_sSq, double &dU_ov_dm, double &dU_ov_ds);
	
	
	//computes gradient and function
	void computeConstrainedEffectivePotential_FunctionAndGradient(const double magnetization, const double staggeredMagnetization, double &U, double &dU_ov_dm, double &dU_ov_ds);
	//versions of computing function and gradient
	void computeFermionicContribution_FunctionAndGradient_qad(const double magnetization, const double staggeredMagnetization, double &Uf, double &dUf_ov_dm, double &dUf_ov_ds);
	//faster
	void computeFermionicContribution_FunctionAndGradient_FromStoredEigenvalues( const double magnetization, const double staggeredMagnetization, double &Uf, double &dUf_ov_dm, double &dUf_ov_ds);
	// format gsl likes
	void computeConstrainedEffectivePotential_FunctionAndGradient_gsl(const gsl_vector *mags, double *U, gsl_vector *gradient_of_U);
	//for nicer code
	void fermionicContributionInline_FunctionAndGradient_FromStoredEigenvalues(int index, double ySq_mSq, double ySq_sSq, double &Uf, double &dUf_ov_dm, double &dUf_ov_ds);
	
	//the second derivative
	void computeConstrainedEffectivePotential_secondDerivatives(const double magnetization, const double staggeredMagnetization, double &d2U_ov_dmdm, double &d2U_ov_dsds, double &d2U_ov_dmds);
	
	
	void compute_fermionicContribution_secondDerivatives_FromStoredEigenvalues(const double magnetization, const double staggeredMagnetization, double &d2Uf_ov_dmdm, double &d2Uf_ov_dsds, double &d2Uf_ov_dmds);
	
	void fermionicContributionInline_secondDerivatives_FromStoredEigenvalues(int index, const double ySq_mSq, const double ySq_sSq, double &d2Uf_ov_dmdm, double &d2Uf_ov_dsds, double &d2Uf_ov_dmds);

	public: 
	
	int initializeMinimizer(const double magnetization, const double staggeredMagnetization);
	
	//will basically do the same as initiallizeMinimizer. Only works if it is already initiallized
	//will be called, if one of the parameters is reset
	int reInitializeMinimizer();//takes current status
	int reInitializeMinimizer(const double magnetization, const double staggeredMagnetization);
	
	
	int iterateMinimizer();
	int itarateUntilToleranceReached(); //max number of iteration given by maxNumerOfIterations
	int itarateUntilToleranceReached(const double tol); //max number of iteration given by maxNumerOfIterations
	int iterateUntilIterationStopps(); //max number of iteration given by maxNumerOfIterations
	
	bool testMinimizerGradient();//will use toleranceForConvergence
	bool testMinimizerGradient(const double tol);
	bool iterationStopped();
	
	
	
	void getActualMinimizerLocation(double &magnetization, double &staggeredMagnetization);
	double getActualMinimizerValue();
	void getActualMinimizerGradient(double &dU_ov_dm, double &dU_ov_ds);
	void getActualSecondDerivative(double &d2U_ov_dmdm, double &d2U_ov_dsds, double &d2U_ov_dmds);
	
	

};

//the wrappers have as a parameter a pointer to the instance of the class
double wrapper_computeConstrainedEffectivePotential_onlyFunction_gsl(const gsl_vector *mags, void *params);
void wrapper_computeConstrainedEffectivePotential_onlyGradient_gsl(const gsl_vector *mags, void *params, gsl_vector *gradient_of_U);
void wrapper_computeConstrainedEffectivePotential_FunctionAndGradient_gsl(const gsl_vector *mags, void *params, double *U, gsl_vector *gradient_of_U);

#endif




