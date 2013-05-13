#include "constrainedEffectivePotential.h"

#include "constrainedEffectivePotential_computeFunctionOnly.cc"
#include "constrainedEffectivePotential_computeGradientOnly.cc"
#include "constrainedEffectivePotential_computeFunctionAndGradient.cc"

constrainedEffectivePotential::constrainedEffectivePotential(const int l0, const int l1, const int l2, const int l3, const bool antiL3):
L0(l0), L1(l1), L2(l2), L3(l3), antiperiodicBC_L3(antiL3),
sinSquaredOfPmu_L0(NULL),sinSquaredOfPmu_L1(NULL),sinSquaredOfPmu_L2(NULL),sinSquaredOfPmu_L3(NULL),
sinSquaredOfPmuHalf_L0(NULL),sinSquaredOfPmuHalf_L1(NULL),sinSquaredOfPmuHalf_L2(NULL),sinSquaredOfPmuHalf_L3(NULL),
cosSquaredOfPmuHalf_L0(NULL),cosSquaredOfPmuHalf_L1(NULL),cosSquaredOfPmuHalf_L2(NULL),cosSquaredOfPmuHalf_L3(NULL),
absNuP(NULL), absNuVarP(NULL), absGammaP(NULL), absGammaVarP(NULL),factorOfMomentum(NULL),
kappa_N(-1.0), lambda_N(-1.0), yukawa_N(-1.0),
N_f(1), rho(1.0), one_ov_twoRho(0.5/rho), r(0.5),
toleranceForLineMinimization(-1.0), toleranceForConvergence(-1.0), initialStepSize(-1.0), maxNumerOfIterations(-1), minimizationAlgorithm(-1), iteratorStoppedFlag(false),
minimizer(NULL)
{
	if( L0%2!=0 || L1%2!=0 || L2%2!=0 || L3%2!=0 )
	{
		std::cerr <<"Error, only even Lattice extends possible!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	fillLatticeMomenta();
	if(L0==L1 && L0==L2){ fillEigenvalues(); }
}


constrainedEffectivePotential::~constrainedEffectivePotential()
{
	if(minimizer!=NULL){ gsl_multimin_fdfminimizer_free (minimizer); }
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
	delete [] absNuP;
	delete [] absNuVarP;
	delete [] absGammaP;
	delete [] absGammaVarP;
	delete [] factorOfMomentum;
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


void constrainedEffectivePotential::fillEigenvalues()
{
	if(!(L0==L1 && L0==L2))
	{
		std::cerr <<"Error, no L^3xL_t lattice in constrainedEffectivePotential::fillEigenvalues()" <<std::endl;
		exit(EXIT_FAILURE);
	}
	int Lhalf=L0/2;
	int LT_half=L3/2;
	int counter=0;
	//do it twice to count
	for(int l3=1; l3<LT_half; ++l3)
	{
		for(int l0=1; l0<Lhalf; ++l0)
		{
			for(int l1=l0+1; l1<Lhalf; ++l1)
			{
				for(int l2=l1+1; l2<Lhalf; ++l2)
				{
					//p,q,r 48
					++counter;
				}
				{
					int l2=l1;
					//p,q,q, 24
					++counter;
					l2=l0;
					//p,p,q, 24
					++counter;
					l2=Lhalf;
					//p,q,L/2 24
					++counter;
					l2=0;
					//0,p,q, 24
					++counter;
				}
			}
			{
				int l1=l0;
				{
					int l2=l1;
					//p,p,p, 8
					++counter;
					l2=Lhalf;
					//p,p,L/2 12
					++counter;
					l2=0;
					//0,p,p, 12
					++counter;
				}
				l1=Lhalf;
				{
					int l2=Lhalf;
					//p,L/2,L/2, 6
					++counter;
					l2=0;
					//0,p,L/2, 12
					++counter;
				}
				l1=0;
				{
					int l2=0;
					//0,0,p  6
					++counter;
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
					++counter;
					l2=0;
					//0,L/2,L/2, 3
					++counter;
				}
				l1=0;
				{
					int l2=0;
					//0,0,L/2, 3
					++counter;
				}
			}
			l0=0;
			{
				int l1=0;
				{
					int l2=0;
					//0,0,0, 1
					++counter;
				}
			}
		}
	}
	//now do the same thing for l3=0 and l3=Lt_half
	for(int l3=0; l3<=LT_half; l3+=(antiperiodicBC_L3 ? LT_half*2 :LT_half))
	{
		for(int l0=1; l0<Lhalf; ++l0)
		{
			for(int l1=l0+1; l1<Lhalf; ++l1)
			{
				for(int l2=l1+1; l2<Lhalf; ++l2)
				{
					//p,q,r 48
					++counter;
				}
				{
					int l2=l1;
					//p,q,q, 24
					++counter;
					l2=l0;
					//p,p,q, 24
					++counter;
					l2=Lhalf;
					//p,q,L/2 24
					++counter;
					l2=0;
					//0,p,q, 24
					++counter;
				}
			}
			{
				int l1=l0;
				{
					int l2=l1;
					//p,p,p, 8
					++counter;
					l2=Lhalf;
					//p,p,L/2 12
					++counter;
					l2=0;
					//0,p,p, 12
					++counter;
				}
				l1=Lhalf;
				{
					int l2=Lhalf;
					//p,L/2,L/2, 6
					++counter;
					l2=0;
					//0,p,L/2, 12
					++counter;
				}
				l1=0;
				{
					int l2=0;
					//0,0,p  6
					++counter;
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
					++counter;
					l2=0;
					//0,L/2,L/2, 3
					++counter;
				}
				l1=0;
				{
					int l2=0;
					//0,0,L/2, 3
					++counter;
				}
			}
			l0=0;
			{
				int l1=0;
				{
					int l2=0;
					//0,0,0, 1
					++counter;
				}
			}
		}
	}
	numberOfDistingtMomenta=counter;
	std::cout <<"numberOfDistingtMomenta: "<<numberOfDistingtMomenta <<std::endl;
	absNuP = new double [numberOfDistingtMomenta];
	absNuVarP = new double [numberOfDistingtMomenta];
	absGammaP = new double [numberOfDistingtMomenta];
	absGammaVarP = new double [numberOfDistingtMomenta];
	factorOfMomentum = new double [numberOfDistingtMomenta];
	std::complex< double > nuOfP(0.,0.), nuOfVarP(0.,0.);
	//factor has to be doubled due to l3
	counter=0;
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
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 96.0;
					++counter;
				}
				{
					int l2=l1;
					//p,q,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 48.0;
					++counter;
					l2=l0;
					//p,p,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 48.0;
					++counter;
					l2=Lhalf;
					//p,q,L/2 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 48.0;
					++counter;
					l2=0;
					//0,p,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 48.0;
					++counter;
				}
			}
			{
				int l1=l0;
				{
					int l2=l1;
					//p,p,p, 8
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 16.0;
					++counter;
					l2=Lhalf;
					//p,p,L/2 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 24.0;
					++counter;
					l2=0;
					//0,p,p, 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 24.0;
					++counter;
				}
				l1=Lhalf;
				{
					int l2=Lhalf;
					//p,L/2,L/2, 6
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 12.0;
					++counter;
					l2=0;
					//0,p,L/2, 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 24.0;
					++counter;
				}
				l1=0;
				{
					int l2=0;
					//0,0,p  6
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 12.0;
					++counter;
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
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 2.0;
					++counter;
					l2=0;
					//0,L/2,L/2, 3
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 6.0;
					++counter;
				}
				l1=0;
				{
					int l2=0;
					//0,0,L/2, 3
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 6.0;
					++counter;
				}
			}
			l0=0;
			{
				int l1=0;
				{
					int l2=0;
					//0,0,0, 1
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = 2.0;
					++counter;
				}
			}
		}
	}
	double BC_factor=(antiperiodicBC_L3 ? 2.0 :1.0);
	//now do the same thing for l3=0 and l3=Lt_half
	for(int l3=0; l3<=LT_half; l3+=(antiperiodicBC_L3 ? LT_half*2.0 :LT_half))
	{
		for(int l0=1; l0<Lhalf; ++l0)
		{
			for(int l1=l0+1; l1<Lhalf; ++l1)
			{
				for(int l2=l1+1; l2<Lhalf; ++l2)
				{
					//p,q,r 48
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*48.0;
					++counter;
				}
				{
					int l2=l1;
					//p,q,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*24.0;
					++counter;
					l2=l0;
					//p,p,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*24.0;
					++counter;
					l2=Lhalf;
					//p,q,L/2 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*24.0;
					++counter;
					l2=0;
					//0,p,q, 24
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*24.0;
					++counter;
				}
			}
			{
				int l1=l0;
				{
					int l2=l1;
					//p,p,p, 8
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*8.0;
					++counter;
					l2=Lhalf;
					//p,p,L/2 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*12.0;
					++counter;
					l2=0;
					//0,p,p, 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*12.0;
					++counter;
				}
				l1=Lhalf;
				{
					int l2=Lhalf;
					//p,L/2,L/2, 6
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*6.0;
					++counter;
					l2=0;
					//0,p,L/2, 12
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*12.0;
					++counter;
				}
				l1=0;
				{
					int l2=0;
					//0,0,p  6
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*6.0;
					++counter;
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
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*1.0;
					++counter;
					l2=0;
					//0,L/2,L/2, 3
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*3.0;
					++counter;
				}
				l1=0;
				{
					int l2=0;
					//0,0,L/2, 3
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*3.0;
					++counter;
				}
			}
			l0=0;
			{
				int l1=0;
				{
					int l2=0;
					//0,0,0, 1
					computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
					absNuP[counter]=std::abs(nuOfP); 
					absNuVarP[counter]=std::abs(nuOfVarP);
					absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
					absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
					factorOfMomentum[counter] = BC_factor*1.0;
					++counter;
				}
			}
		}
	}
}

void constrainedEffectivePotential::set_kappa_N(double new_k){ kappa_N=new_k; }
void constrainedEffectivePotential::set_lambda_N(double new_l){ lambda_N=new_l; }
void constrainedEffectivePotential::set_yukawa_N(double new_y){ yukawa_N=new_y; }

void constrainedEffectivePotential::set_N_f(int new_N){ N_f=new_N; }
void constrainedEffectivePotential::set_rho(double new_rho){ rho=new_rho; one_ov_twoRho=0.5/rho; }
void constrainedEffectivePotential::set_r(double new_r){ r=new_r; }

void constrainedEffectivePotential::set_toleranceForLineMinimization(double new_tol){ toleranceForLineMinimization=new_tol; }
void constrainedEffectivePotential::set_toleranceForConvergence(double new_tol){ toleranceForConvergence=new_tol; }
void constrainedEffectivePotential::set_initialStepSize(double new_step){ initialStepSize=new_step; }
void constrainedEffectivePotential::set_maxNumerOfIterations(int new_NOI){ maxNumerOfIterations=new_NOI; }
void constrainedEffectivePotential::set_minimizationAlgorithm(int new_alg)
{
// 	delete algorithmForMinimization;
// 	algorithmForMinimization = new gsl_multimin_fdfminimizer_type;
	switch(new_alg)
	{
		case 1: 
			minimizationAlgorithm=new_alg;
			algorithmForMinimization = (gsl_multimin_fdfminimizer_conjugate_fr);
			break;
		case 2: 
			minimizationAlgorithm=new_alg;
			algorithmForMinimization = (gsl_multimin_fdfminimizer_conjugate_pr);
			break;
		case 3: 
			minimizationAlgorithm=new_alg;
			algorithmForMinimization = (gsl_multimin_fdfminimizer_vector_bfgs2);
			break;
		case 4: 
			minimizationAlgorithm=new_alg;
			algorithmForMinimization = (gsl_multimin_fdfminimizer_vector_bfgs);
			break;
		case 5: 
			minimizationAlgorithm=new_alg;
			algorithmForMinimization = (gsl_multimin_fdfminimizer_steepest_descent);
			break;
		default:
			std::cerr <<"Error, algorithm " <<new_alg <<" is not valid in constrainedEffectivePotential::set_minimizationAlgorithm(int new_alg) " <<std::endl;
			std::cerr <<"Use one of the following:" <<std::endl;
			std::cerr <<"1 = gsl_multimin_fdfminimizer_conjugate_fr" <<std::endl;
			std::cerr <<"2 = gsl_multimin_fdfminimizer_conjugate_pr" <<std::endl;
			std::cerr <<"3 = gsl_multimin_fdfminimizer_vector_bfgs2" <<std::endl;
			std::cerr <<"4 = gsl_multimin_fdfminimizer_vector_bfgs" <<std::endl;
			std::cerr <<"5 = gsl_multimin_fdfminimizer_steepest_descent" <<std::endl;
			exit(EXIT_FAILURE);
	}
}


std::complex< double > constrainedEffectivePotential::computeAnalyticalEigenvalue(const double p0, const double p1, const double p2, const double p3)
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
	
	double one_ov_denom=1.0/sqrt(p_tilde_sq + (r*p_hat_sq - rho)*(r*p_hat_sq - rho) );
	
	return std::complex< double >(rho + rho*one_ov_denom*(r*p_hat_sq - rho) , rho*one_ov_denom*sqrt(p_tilde_sq) );
}



std::complex< double > constrainedEffectivePotential::computeAnalyticalEigenvalue_fromIndex(const int l0, const int l1, const int l2, const int l3)
{
	double p_tilde_sq = sinSquaredOfPmu_L0[l0] + sinSquaredOfPmu_L1[l1] + 
	                    sinSquaredOfPmu_L2[l2] + sinSquaredOfPmu_L3[l3];
	
	
	double p_hat_sq = sinSquaredOfPmuHalf_L0[l0] + sinSquaredOfPmuHalf_L1[l1] +
	                  sinSquaredOfPmuHalf_L2[l2] + sinSquaredOfPmuHalf_L3[l3];
	p_hat_sq*=4.0;

	
	double one_ov_denom=1.0/sqrt(p_tilde_sq + (r*p_hat_sq - rho)*(r*p_hat_sq - rho) );
	return std::complex< double >(rho + rho*one_ov_denom*(r*p_hat_sq - rho) , rho*one_ov_denom*sqrt(p_tilde_sq) );
}



void constrainedEffectivePotential::computeAnalyticalEigenvalue_fromIndex_pAndVarP(const int l0, const int l1, const int l2, const int l3, std::complex< double > &nuOfP, std::complex< double > &nuOfVarP)
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



int constrainedEffectivePotential::initializeMinimizer(const double magnetization, const double staggeredMagnetization)
{
	if( toleranceForLineMinimization < 0.0)
	{
		std::cerr <<"Error, no or negative tolerance set in constrainedEffectivePotential::initializeMinimizer(...)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if( initialStepSize <= 0.0)
	{
		std::cerr <<"Error, no or negative tolerance set in constrainedEffectivePotential::initializeMinimizer(...)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	if(minimizationAlgorithm<1 || minimizationAlgorithm >5 || algorithmForMinimization==NULL)
	{
		std::cerr <<"No valid minimizationAlgorithm algorithm set in constrainedEffectivePotential::initializeMinimizer(...)" <<std::endl;
		std::cerr <<"Use one of the following:" <<std::endl;
		std::cerr <<"1 = gsl_multimin_fdfminimizer_conjugate_fr" <<std::endl;
		std::cerr <<"2 = gsl_multimin_fdfminimizer_conjugate_pr" <<std::endl;
		std::cerr <<"3 = gsl_multimin_fdfminimizer_vector_bfgs2" <<std::endl;
		std::cerr <<"4 = gsl_multimin_fdfminimizer_vector_bfgs" <<std::endl;
		std::cerr <<"5 = gsl_multimin_fdfminimizer_steepest_descent" <<std::endl;
		exit(EXIT_FAILURE);
	}
	
	gsl_vector *mags = gsl_vector_alloc(2);
	gsl_vector_set(mags, 0, magnetization);
	gsl_vector_set(mags, 1, staggeredMagnetization);
	//allocate minimizer
	minimizer = gsl_multimin_fdfminimizer_alloc (algorithmForMinimization, 2);
	//set the handler
	functionHandler.n=2;
	functionHandler.f = &wrapper_computeConstrainedEffectivePotential_onlyFunction_gsl;
	functionHandler.df= &wrapper_computeConstrainedEffectivePotential_onlyGradient_gsl;
	functionHandler.fdf=&wrapper_computeConstrainedEffectivePotential_FunctionAndGradient_gsl;
	functionHandler.params=(void *) this;
	
	int ret = gsl_multimin_fdfminimizer_set(minimizer, &functionHandler, mags, initialStepSize, toleranceForLineMinimization);
	gsl_vector_free(mags);
	return ret;
}


int constrainedEffectivePotential::iterateMinimizer()
{
	iteratorStoppedFlag=false;
	int returnValue=gsl_multimin_fdfminimizer_iterate(minimizer);
	if(returnValue==GSL_ENOPROG){iteratorStoppedFlag=true;}
	return returnValue;
}

bool constrainedEffectivePotential::testMinimizerGradient()//will use toleranceForConvergence
{
	return testMinimizerGradient(toleranceForConvergence);
}

bool constrainedEffectivePotential::testMinimizerGradient(const double tol)
{
	int testResult=gsl_multimin_test_gradient( gsl_multimin_fdfminimizer_gradient (minimizer), tol);
	if(testResult==GSL_SUCCESS){return true;}
	else{return false;}
}

bool constrainedEffectivePotential::iterationStopped()
{
	return iteratorStoppedFlag;
}

void constrainedEffectivePotential::getActualMinimizerLocation(double &magnetization, double &staggeredMagnetization)
{
	magnetization=gsl_vector_get( gsl_multimin_fdfminimizer_x(minimizer),0);
	staggeredMagnetization=gsl_vector_get( gsl_multimin_fdfminimizer_x(minimizer),1);	
}

double constrainedEffectivePotential::getActualMinimizerValue()
{
	return gsl_multimin_fdfminimizer_minimum(minimizer);
}

void constrainedEffectivePotential::getActualMinimizerGradient(double &dU_ov_dm, double &dU_ov_ds)
{
	dU_ov_dm=gsl_vector_get( gsl_multimin_fdfminimizer_gradient(minimizer),0);
	dU_ov_ds=gsl_vector_get( gsl_multimin_fdfminimizer_gradient(minimizer),1);	
}






