#include "constrainedEffectivePotential.h"

#include "constrainedEffectivePotential_computeFunctionOnly.cc"
#include "constrainedEffectivePotential_computeGradientOnly.cc"
#include "constrainedEffectivePotential_computeFunctionAndGradient.cc"
#include "constrainedEffectivePotential_computeBosonicLoop.cc"
#include "constrainedEffectivePotential_computeSecondDerivative.cc"

constrainedEffectivePotential::constrainedEffectivePotential(const int l0, const int l1, const int l2, const int l3, const bool antiL3):
L0(l0), L1(l1), L2(l2), L3(l3), antiperiodicBC_L3(antiL3),
sinSquaredOfPmu_L0(NULL),sinSquaredOfPmu_L1(NULL),sinSquaredOfPmu_L2(NULL),sinSquaredOfPmu_L3(NULL),
sinSquaredOfPmuHalf_L0(NULL),sinSquaredOfPmuHalf_L1(NULL),sinSquaredOfPmuHalf_L2(NULL),sinSquaredOfPmuHalf_L3(NULL),
cosSquaredOfPmuHalf_L0(NULL),cosSquaredOfPmuHalf_L1(NULL),cosSquaredOfPmuHalf_L2(NULL),cosSquaredOfPmuHalf_L3(NULL),
cosOfPmu_L0(NULL), cosOfPmu_L1(NULL), cosOfPmu_L2(NULL), cosOfPmu_L3(NULL),
numberOfDistingtMomenta(-1), absNuP(NULL), absNuVarP(NULL), absGammaP(NULL), absGammaVarP(NULL),factorOfMomentum(NULL),
numberOfDistingtMomenta_bosonic(-1), sumOfCosOfPmu(NULL), factorOfMomentum_bosonic(NULL),
kappa_N(-1.0), lambda_N(-1.0), yukawa_N(-1.0),
N_f(1), rho(1.0), one_ov_twoRho(0.5/rho), r(0.5),
bosonicLoop(-1.0), bosonicLoopSet(false), useBosonicLoop(true),
toleranceForLineMinimization(-1.0), toleranceForConvergence(-1.0), initialStepSize(-1.0), maxNumerOfIterations(100), minimizationAlgorithm(-1), minimizerInitialized(false), iteratorStoppedFlag(false),
minimizer(NULL)
{
	if( L0%2!=0 || L1%2!=0 || L2%2!=0 || L3%2!=0 )
	{
		std::cerr <<"Error, only even Lattice extends possible!" <<std::endl;
		exit(EXIT_FAILURE);
	}
	fillLatticeMomenta();
	fillEigenvalues();
	fillBosonicLoopValues();
}


constrainedEffectivePotential::~constrainedEffectivePotential()
{
	if(minimizer!=NULL){ gsl_multimin_fdfminimizer_free (minimizer); }
	delete [] sinSquaredOfPmu_L0; delete [] sinSquaredOfPmu_L1; 
	delete [] sinSquaredOfPmu_L2; delete [] sinSquaredOfPmu_L3;
	delete [] sinSquaredOfPmuHalf_L0; delete [] sinSquaredOfPmuHalf_L1;
	delete [] sinSquaredOfPmuHalf_L2; delete [] sinSquaredOfPmuHalf_L3;
	delete [] cosSquaredOfPmuHalf_L0; delete [] cosSquaredOfPmuHalf_L1;
	delete [] cosSquaredOfPmuHalf_L2; delete [] cosSquaredOfPmuHalf_L3;
	delete [] cosOfPmu_L0; delete [] cosOfPmu_L1;
	delete [] cosOfPmu_L2; delete [] cosOfPmu_L3;
	delete [] absNuP; delete [] absNuVarP;
	delete [] absGammaP; delete [] absGammaVarP;
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
	cosOfPmu_L0 = new double [L0]; cosOfPmu_L1 = new double [L1]; cosOfPmu_L2 = new double [L2]; cosOfPmu_L3 = new double [L3];
	double p,sinP,cosP;
	for(int i=0; i<L0; ++i)
	{
		p=2.0 * PI * i * one_ov_L0; 
		sinP=sin(p); sinSquaredOfPmu_L0[i]=sinP*sinP;
		sinP=sin(0.5*p); sinSquaredOfPmuHalf_L0[i]=sinP*sinP;
		cosP=cos(0.5*p); cosSquaredOfPmuHalf_L0[i]=cosP*cosP;
		cosP=cos(p); cosOfPmu_L0[i]=cosP;
// 		std::cout <<"l0 = " << i <<"   p_mu = " <<p <<"   sinSquaredOfPmu = " <<sinSquaredOfPmu_L0[i] <<"   sinSquaredOfPmuHalf = " <<sinSquaredOfPmuHalf_L0[i] <<std::endl;
	}
	for(int i=0; i<L1; ++i)
	{
		p=2.0 * PI * i * one_ov_L1; 
		sinP=sin(p); sinSquaredOfPmu_L1[i]=sinP*sinP;
		sinP=sin(0.5*p); sinSquaredOfPmuHalf_L1[i]=sinP*sinP;
		cosP=cos(0.5*p); cosSquaredOfPmuHalf_L1[i]=cosP*cosP;
		cosP=cos(p); cosOfPmu_L1[i]=cosP;
// 		std::cout <<"l1 = " << i <<"   p_mu = " <<p <<"   sinSquaredOfPmu = " <<sinSquaredOfPmu_L1[i] <<"   sinSquaredOfPmuHalf = " <<sinSquaredOfPmuHalf_L1[i] <<std::endl;
	}
	for(int i=0; i<L2; ++i)
	{
		p=2.0 * PI * i * one_ov_L2; 
		sinP=sin(p); sinSquaredOfPmu_L2[i]=sinP*sinP;
		sinP=sin(0.5*p); sinSquaredOfPmuHalf_L2[i]=sinP*sinP;
		cosP=cos(0.5*p); cosSquaredOfPmuHalf_L2[i]=cosP*cosP;
		cosP=cos(p); cosOfPmu_L2[i]=cosP;
// 		std::cout <<"l2 = " << i <<"   p_mu = " <<p <<"   sinSquaredOfPmu = " <<sinSquaredOfPmu_L2[i] <<"   sinSquaredOfPmuHalf = " <<sinSquaredOfPmuHalf_L2[i] <<std::endl;
	}
	for(int i=0; i<L3; ++i)
	{
		p=2.0 * PI * i * one_ov_L3 + toAddForL3; 
		sinP=sin(p); sinSquaredOfPmu_L3[i]=sinP*sinP;
		sinP=sin(0.5*p); sinSquaredOfPmuHalf_L3[i]=sinP*sinP;
		cosP=cos(0.5*p); cosSquaredOfPmuHalf_L3[i]=cosP*cosP;
		p=2.0 * PI * i * one_ov_L3;
		cosP=cos(p); cosOfPmu_L3[i]=cosP;
// 		std::cout <<"l3 = " << i <<"   p_mu = " <<p <<"   sinSquaredOfPmu = " <<sinSquaredOfPmu_L3[i] <<"   sinSquaredOfPmuHalf = " <<sinSquaredOfPmuHalf_L3[i] <<std::endl;
	}
}


void constrainedEffectivePotential::fillEigenvalues()
{
	delete [] absNuP;
	delete [] absNuVarP;
	delete [] absGammaP;
	delete [] absGammaVarP;
	delete [] factorOfMomentum;
	if(!(L0==L1 && L0==L2))
	{
		int L0_half=L0/2;
		int L1_half=L1/2;
		int L2_half=L2/2;
		int L3_half=L3/2;
		int dummyVolume(0);
		if(antiperiodicBC_L3)
		{
			dummyVolume=(L0_half+1)*(L1_half+1)*(L2_half+1)*(L3_half);
		}
		else
		{
			dummyVolume=(L0_half+1)*(L1_half+1)*(L2_half+1)*(L3_half+1);
		}
		numberOfDistingtMomenta=dummyVolume;
		std::cout <<"numberOfDistingtMomenta: "<<numberOfDistingtMomenta <<std::endl;
		absNuP = new double [dummyVolume];
		absNuVarP = new double [dummyVolume];
		absGammaP = new double [dummyVolume];
		absGammaVarP = new double [dummyVolume];
		factorOfMomentum = new double [dummyVolume];
		double fac_l0,fac_l1,fac_l2,fac_l3;
		std::complex< double > nuOfP(0.,0.), nuOfVarP(0.,0.);
		size_t counter=0;
		for(int l0=0; l0<=L0_half; ++l0)
		{
			if(l0==0 || l0==L0_half){ fac_l0=1.0; }else{ fac_l0=2.0; }
			for(int l1=0; l1<=L1_half; ++l1)
			{
				if(l1==0 || l1==L1_half){ fac_l1=1.0; }else{ fac_l1=2.0; }
				for(int l2=0; l2<=L2_half; ++l2)
				{
					if(l2==0 || l2==L2_half){ fac_l2=1.0; }else{ fac_l2=2.0; }
					for(int l3=0; l3<=L3_half; ++l3)
					{
						if(!antiperiodicBC_L3 &&( l3==0 || l3==L3_half) ){ fac_l3=1.0; }
						else if(!antiperiodicBC_L3 && !( l3==0 || l3==L3_half) ){ fac_l3=2.0; }
						else if(antiperiodicBC_L3 && l3==L3_half){ continue; }
						else{ fac_l3=2.0; }
						
						computeAnalyticalEigenvalue_fromIndex_pAndVarP( l0,l1,l2,l3, nuOfP, nuOfVarP );
						absNuP[counter]=std::abs(nuOfP); 
						absNuVarP[counter]=std::abs(nuOfVarP);
						absGammaP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfP );
						absGammaVarP[counter]=std::abs( 1.0 - one_ov_twoRho*nuOfVarP );
						factorOfMomentum[counter] = fac_l0*fac_l1*fac_l2*fac_l3;
						++counter;
					}
				}
			}
		}
	}
	else
	{
		int Lhalf=L0/2;
		int LT_half=L3/2;
		size_t counter=0;
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
						l2=l2; //no compiler warning
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
						l2=l2; //no compiler warning
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
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2; //no compiler warning
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
						l2=l2; //no compiler warning
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
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2;//no compiler warning
						//0,0,0, 1
						++counter;
					}
				}
			}
		}
		numberOfDistingtMomenta=counter;
		std::cout <<"numberOfDistingtMomenta: "<<numberOfDistingtMomenta <<std::endl;
		//delete first
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
}




void constrainedEffectivePotential::fillBosonicLoopValues()
{
	delete [] sumOfCosOfPmu; delete [] factorOfMomentum_bosonic;
	if(  !(L0==L1 && L0==L2))
	{
		int L0_half=L0/2, L1_half=L1/2, L2_half=L2/2, L3_half=L3/2;
		int dummyVolume=(L0_half+1)*(L1_half+1)*(L2_half+1)*(L3_half+1)-2;
		numberOfDistingtMomenta_bosonic=dummyVolume;
		std::cout <<"numberOfDistingtMomenta_bosonic: "<<numberOfDistingtMomenta_bosonic <<std::endl;
		sumOfCosOfPmu = new double [dummyVolume];
		factorOfMomentum_bosonic = new double [dummyVolume];
		double fac_l0,fac_l1,fac_l2,fac_l3;
		size_t counter=0;
		for(int l0=0; l0<=L0_half; ++l0)
		{
			if(l0==0 || l0==L0_half){ fac_l0=1.0; }else{ fac_l0=2.0; }
			for(int l1=0; l1<=L1_half; ++l1)
			{
				if(l1==0 || l1==L1_half){ fac_l1=1.0; }else{ fac_l1=2.0; }
				for(int l2=0; l2<=L2_half; ++l2)
				{
					if(l2==0 || l2==L2_half){ fac_l2=1.0; }else{ fac_l2=2.0; }
					for(int l3=((l0+l1+l2)?(0):(1)); l3<=( (l0==L0_half && l1==L1_half && l2==L2_half)?(L3_half-1):L3_half ); ++l3)
					{
						if(l3==0 || l3==L3_half){ fac_l3=1.0; }else{ fac_l3=2.0; }
						
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = fac_l0*fac_l1*fac_l2*fac_l3;
						++counter;
					}
				}
			}
		}
		std::cout <<"Momenta entries in sumOfCosOfPmu: " <<counter <<std::endl;
	}
	else
	{
		int Lhalf=L0/2;
		int LT_half=L3/2;
		size_t counter=0;
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
						l2=l2; //no compiler warning
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
						l2=l2; //no compiler warning
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
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,0, 1
						++counter;
					}
				}
			}
		}
		//now do the same thing for l3=0 and l3=Lt_half
		for(int l3=0; l3<=LT_half; l3+=LT_half )
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
						l2=l2; //no compiler warning
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
						if(l3!=LT_half){++counter; }
						l2=0;
						//0,L/2,L/2, 3
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2;//no compiler warning
						//0,0,0, 1
						if(l3 != 0){++counter;}
					}
				}
			}
		}
		numberOfDistingtMomenta_bosonic=counter;
		std::cout <<"numberOfDistingtMomenta_bosonic: "<<numberOfDistingtMomenta_bosonic <<std::endl;
		sumOfCosOfPmu = new double [numberOfDistingtMomenta_bosonic];
		factorOfMomentum_bosonic = new double [numberOfDistingtMomenta_bosonic];
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
						l2=l2; //no compiler warning
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 96.0;
						++counter;
					}
					{
						int l2=l1;
						//p,q,q, 24
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 48.0;
						++counter;
						l2=l0;
						//p,p,q, 24
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 48.0;
						++counter;
						l2=Lhalf;
						//p,q,L/2 24
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 48.0;
						++counter;
						l2=0;
						//0,p,q, 24
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 48.0;
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						//p,p,p, 8
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 16.0;
						++counter;
						l2=Lhalf;
						//p,p,L/2 12
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 24.0;
						++counter;
						l2=0;
						//0,p,p, 12
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 24.0;
						++counter;
					}
					l1=Lhalf;
					{
						int l2=Lhalf;
						//p,L/2,L/2, 6
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 12.0;
						++counter;
						l2=0;
						//0,p,L/2, 12
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 24.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,p  6
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 12.0;
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
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 2.0;
						++counter;
						l2=0;
						//0,L/2,L/2, 3
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 6.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 6.0;
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,0, 1
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 2.0;
						++counter;
					}
				}
			}
		}
		//now do the same thing for l3=0 and l3=Lt_half
		for(int l3=0; l3<=LT_half; l3+=LT_half )
		{
			for(int l0=1; l0<Lhalf; ++l0)
			{
				for(int l1=l0+1; l1<Lhalf; ++l1)
				{
					for(int l2=l1+1; l2<Lhalf; ++l2)
					{
						//p,q,r 48
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 48.0;
						++counter;
					}
					{
						int l2=l1;
						//p,q,q, 24
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 24.0;
						++counter;
						l2=l0;
						//p,p,q, 24
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 24.0;
						++counter;
						l2=Lhalf;
						//p,q,L/2 24
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 24.0;
						++counter;
						l2=0;
						//0,p,q, 24
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 24.0;
						++counter;
					}
				}
				{
					int l1=l0;
					{
						int l2=l1;
						//p,p,p, 8
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 8.0;
						++counter;
						l2=Lhalf;
						//p,p,L/2 12
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 12.0;
						++counter;
						l2=0;
						//0,p,p, 12
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 12.0;
						++counter;
					}
					l1=Lhalf;
					{
						int l2=Lhalf;
						//p,L/2,L/2, 6
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 6.0;
						++counter;
						l2=0;
						//0,p,L/2, 12
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 12.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,p  6
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 6.0;
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
						if(l3!=LT_half)
						{
							sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
							factorOfMomentum_bosonic[counter] = 1.0;
							++counter; 
						}
						l2=0;
						//0,L/2,L/2, 3
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 3.0;
						++counter;
					}
					l1=0;
					{
						int l2=0;
						l2=l2; //no compiler warning
						//0,0,L/2, 3
						sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
						factorOfMomentum_bosonic[counter] = 3.0;
						++counter;
					}
				}
				l0=0;
				{
					int l1=0;
					l1=l1; //no compiler warning
					{
						int l2=0;
						l2=l2;//no compiler warning
						//0,0,0, 1
						if(l3 != 0)
						{
							sumOfCosOfPmu[counter]=cosOfPmu_L0[l0]+cosOfPmu_L1[l1]+cosOfPmu_L2[l2]+cosOfPmu_L3[l3];
							factorOfMomentum_bosonic[counter] = 1.0;
							++counter;
						}
					}
				}
			}
		}
	}
		
	

}

void constrainedEffectivePotential::set_kappa_N(double new_k)
{
	if( new_k==kappa_N ){ return; }
	kappa_N=new_k; reInitializeMinimizer();
	if(bosonicLoopSet && useBosonicLoop)
	{
		bosonicLoop=computeBosonicPropagatorSum_fromStoredSumOfCos();
	}
}
void constrainedEffectivePotential::set_lambda_N(double new_l)
{ 
	if( new_l==lambda_N ){ return; }
	lambda_N=new_l; reInitializeMinimizer();
	if(bosonicLoopSet && useBosonicLoop)
	{
		bosonicLoop=computeBosonicPropagatorSum_fromStoredSumOfCos();
	}
}
void constrainedEffectivePotential::set_yukawa_N(double new_y){ yukawa_N=new_y; reInitializeMinimizer();}
void constrainedEffectivePotential::set_kappa_lambda_yukawa_N(double new_k, double new_l, double new_y)
{
	bool change=(new_l!=lambda_N || new_k!=kappa_N);
	kappa_N=new_k; lambda_N=new_l; yukawa_N=new_y; reInitializeMinimizer();
	if(bosonicLoopSet && useBosonicLoop && change)
	{
		bosonicLoop=computeBosonicPropagatorSum_fromStoredSumOfCos();
	}
}

void constrainedEffectivePotential::set_N_f(int new_N){ N_f=new_N; }
void constrainedEffectivePotential::set_rho(double new_rho){ rho=new_rho; one_ov_twoRho=0.5/rho; fillEigenvalues(); reInitializeMinimizer(); }
void constrainedEffectivePotential::set_r(double new_r){ r=new_r; fillEigenvalues(); reInitializeMinimizer(); }
void constrainedEffectivePotential::set_useBosonicLoop(bool newSet){ useBosonicLoop=newSet; }

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


double constrainedEffectivePotential::get_kappa_N(){ return kappa_N; }
double constrainedEffectivePotential::get_lambda_N(){ return lambda_N; }
double constrainedEffectivePotential::get_yukawa_N(){ return yukawa_N; }

int constrainedEffectivePotential::get_N_f(){ return N_f; }
double constrainedEffectivePotential::get_rho(){ return rho; }
double constrainedEffectivePotential::get_r(){ return r; }

double  constrainedEffectivePotential::get_toleranceForLineMinimization(){ return toleranceForLineMinimization; }
double  constrainedEffectivePotential::get_toleranceForConvergence(){ return toleranceForConvergence; }
double  constrainedEffectivePotential::get_initialStepSize(){ return initialStepSize; }
int  constrainedEffectivePotential::get_maxNumerOfIterations(){ return maxNumerOfIterations; }
int  constrainedEffectivePotential::get_minimizationAlgorithm(){ return minimizationAlgorithm; }



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
	if(minimizerInitialized){ return reInitializeMinimizer(magnetization, staggeredMagnetization); }
	
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
	minimizerInitialized=true;
	return ret;
}



int constrainedEffectivePotential::reInitializeMinimizer(const double magnetization, const double staggeredMagnetization)
{
	if(!minimizerInitialized){ return 1; }
	gsl_multimin_fdfminimizer_free( minimizer );
	minimizer = gsl_multimin_fdfminimizer_alloc (algorithmForMinimization, 2);
	gsl_vector *mags = gsl_vector_alloc(2);
	gsl_vector_set(mags, 0, magnetization);
	gsl_vector_set(mags, 1, staggeredMagnetization);
	int ret = gsl_multimin_fdfminimizer_set(minimizer, &functionHandler, mags, initialStepSize, toleranceForLineMinimization);
	gsl_vector_free(mags);
	iteratorStoppedFlag=false;
	return ret;
}



int constrainedEffectivePotential::reInitializeMinimizer()
{
	if(!minimizerInitialized){ return 1; }
	double mag(0.0),stag(0.0);
	getActualMinimizerLocation(mag,stag);
	return reInitializeMinimizer(mag, stag);
}




int constrainedEffectivePotential::iterateMinimizer()
{
	iteratorStoppedFlag=false;
	int returnValue=gsl_multimin_fdfminimizer_iterate(minimizer);
	if(returnValue==GSL_ENOPROG){iteratorStoppedFlag=true;}
	return returnValue;
}

int constrainedEffectivePotential::itarateUntilToleranceReached()
{
	return itarateUntilToleranceReached(toleranceForConvergence);
}

int constrainedEffectivePotential::itarateUntilToleranceReached(const double tol)
{
	if(!minimizerInitialized)
	{
		std::cerr <<"Error, minimizer not initiallized in constrainedEffectivePotential::itarateUntilToleranceReached(const double tol)" <<std::endl;
		exit(EXIT_FAILURE);
	}
	
	for(int i=1; i<=maxNumerOfIterations; ++i)
	{
		iterateMinimizer();
		if(testMinimizerGradient(tol)){ return i; }
		if( iterationStopped() ){ return i; }
	}
	return 0;
}


int constrainedEffectivePotential::iterateUntilIterationStopps()
{
	if(!minimizerInitialized)
	{
		std::cerr <<"Error, minimizer not initiallized in constrainedEffectivePotential::iterateUntilIterationStopps()" <<std::endl;
		exit(EXIT_FAILURE);
	}
	
	for(int i=1; i<=maxNumerOfIterations; ++i)
	{
		iterateMinimizer();
		if( iterationStopped() ){ return i; }
	}
	return 0;
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

void constrainedEffectivePotential::getActualSecondDerivative(double &d2U_ov_dmdm, double &d2U_ov_dsds, double &d2U_ov_dmds)
{
	double m,s;
	getActualMinimizerLocation(m, s);
	computeConstrainedEffectivePotential_secondDerivatives(m, s, d2U_ov_dmdm, d2U_ov_dsds, d2U_ov_dmds);
}




