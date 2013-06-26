#include "constrainedEffectivePotential.h"

/*
double constrainedEffectivePotential::computeBosonicPropagatorSum_qad()
{
	double twoMinFourLambdaN=2.0-4.0*lambda_N;
	double fourKappa=4.0*kappa_N;
	double dummyForAddition(0.0);
// 	double counter=0;
	int L0_half=L0/2, L1_half=L1/2, L2_half=L2/2, L3_half=L3/2;
	for(int l0=0; l0<L0; ++l0)
	{
		for(int l1=0; l1<L1; ++l1)
		{
			for(int l2=0; l2<L2; ++l2)
			{
				for(int l3=((l0+l1+l2)?0:1); l3<L3_half; ++l3)
				{
					
					dummyForAddition+=1.0/(twoMinFourLambdaN - fourKappa*(cosOfPmu_L0[l0] + cosOfPmu_L1[l1] + cosOfPmu_L2[l2] + cosOfPmu_L3[l3]));
// 					counter++;
				}
				for(int l3=((l0==L0_half && l1==L1_half && l2==L2_half )?L3_half+1:L3_half); l3<L3; ++l3)
				{
					dummyForAddition+=1.0/(twoMinFourLambdaN - fourKappa*(cosOfPmu_L0[l0] + cosOfPmu_L1[l1] + cosOfPmu_L2[l2] + cosOfPmu_L3[l3]));
// 					counter++;
				}
				
			}
		}
	}
// 	std::cout <<"Counter in constrainedEffectivePotential::computeBosonicPropagatorSum_qad() = " <<counter << std::endl;
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
	return dummyForAddition;
}
*/

double constrainedEffectivePotential::computeBosonicPropagatorSum_fromStoredSumOfCos()
{
	double twoMinFourLambdaN=2.0-4.0*lambda_N;
	double fourKappa=4.0*kappa_N;
	double dummyForAddition(0.0);
	;
	for(int index=0; index<numberOfDistingtMomenta_bosonic; ++index)
	{
		dummyForAddition+=factorOfMomentum_bosonic[index]/(twoMinFourLambdaN - fourKappa*sumOfCosOfPmu[index]);
	}
// 	std::cout <<"Counter in constrainedEffectivePotential::computeBosonicPropagatorSum_qad() = " <<counter << std::endl;
	dummyForAddition/=static_cast< double >(L0); dummyForAddition/=static_cast< double >(L1); dummyForAddition/=static_cast< double >(L2); dummyForAddition/=static_cast< double >(L3); 
	return dummyForAddition;
}
