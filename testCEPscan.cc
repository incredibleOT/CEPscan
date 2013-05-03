#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>

#include "constrainedEffectivePotential.h"


using std::cout;
using std::cerr;
using std::endl;

int main(int narg,char **arg)
{
	
	int L0(4), L1(4), L2(4), L3(6);
// 	int L0(8), L1(10), L2(14), L3(22);
// 	bool antiL3(false);
	bool antiL3(false);
	
	double y=175.0/246.0;
// 	double y=0.0;
	
	double kappa=0.35;
	double lambda=0.01;
	
	constrainedEffectivePotential CEP(L0,L1,L2,L3,antiL3);
	cout <<"constrainedEffectivePotential initialized" <<endl;
	CEP.set_yukawa_N(y);
	CEP.set_kappa_N(kappa);
	CEP.set_lambda_N(lambda);
	
	double mCentral=1.0;
	double sCentral=1.0;
	double mStep=0.0001;
	double sStep=0.0001;
	int steps=5;
	int maxCounter(2*steps+1);
	
	double *resultsMvar = new double [maxCounter];
	double *resultsSvar = new double [maxCounter];
	double m,s;
	
	for(int im=0; im<maxCounter; ++im)
	{
		m=mCentral-mStep*(steps-im);
		resultsMvar[im]=CEP.computeConstrainedEffectivePotential_onlyFunction(m,sCentral);
		s=sCentral-sStep*(steps-im);
		resultsSvar[im]=CEP.computeConstrainedEffectivePotential_onlyFunction(mCentral,s);
	}
	
	for(int i=0; i<maxCounter; ++i)
	{
		cout.precision(15);
		m=mCentral-mStep*(steps-i);
		s=sCentral-sStep*(steps-i);
		cout <<"U(m=" <<m <<", s=" <<sCentral <<") = " <<resultsMvar[i] <<"   U(m=" <<mCentral <<", s=" <<s <<") = " <<resultsSvar[i] <<endl;
	}
	for(int i=0; i<steps; ++i)
	{
		cout.precision(15);
		double dist_m=mStep*(i+1);
		double dist_s=sStep*(i+1);
		cout <<"(U(m=" <<mCentral <<"+" <<dist_m <<", s=" <<sCentral <<") -  U(m=" <<mCentral <<"-" <<dist_m <<", s=" <<sCentral <<"))/"<<2*dist_m <<"  = " 
		     <<(resultsMvar[steps+1+i]-resultsMvar[steps-i-1])/(2.0*dist_m) <<endl;
		     
		cout <<"(U(m=" <<mCentral <<", s=" <<sCentral <<"+" <<dist_s <<") -  U(m=" <<mCentral <<", s=" <<sCentral <<"-" <<dist_s <<"))/"<<2*dist_s <<"  = " 
		     <<(resultsSvar[steps+1+i]-resultsSvar[steps-1-i])/(2.0*dist_s) <<endl;
	}
	
	double res_U_dm(0.0), res_U_ds(0.0);
	CEP.computeConstrainedEffectivePotential_onlyGradient(mCentral,sCentral, res_U_dm, res_U_ds);
	cout <<"dU/dm(m=" <<mCentral <<", s=" <<sCentral <<") = " <<res_U_dm <<"   dU/ds(m=" <<mCentral <<", s=" <<sCentral <<") = " <<res_U_ds <<endl;
	
	delete [] resultsMvar;
	delete [] resultsSvar;
	
	
	


}
