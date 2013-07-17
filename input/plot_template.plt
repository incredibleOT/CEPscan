
L_s                       16
L_t                       32
antiperiodic_L_t          0

#whether or not the first order in lambda should be included(default 1)
include_bosonic_loop       0

#whether an improved treelevel ansatz should be used coming from 
#including contributions from the phi^4 and phi^6 term to the gaussian contribution
use_improved_treeLevel     0

#use first order in lambda for the improved treelevel
use_improved_firstOrder    0

# default: N_f=1, rho=1, r=0.5
# N_f                       1
# rho                       1.0
# r                         0.5

	
# scan_kappa_N            1
kappa_N                   0.100
# kappa_N_min               -0.5
# kappa_N_max               +0.5
# kappa_N_step              +0.01

# scan_lambda_N             1
lambda_N                    0.3
# lambda_N_min                0.0
# lambda_N_max                0.5
# lambda_N_step               0.01

# scan_lambda_6_N             1
lambda_6_N                    0.3
# lambda_6_N_min                0.0
# lambda_6_N_max                0.5
# lambda_6_N_step               0.01

# 175/246=0.711382114    500/246=2.032520325
# scan_yukawa_N           1
yukawa_N                  0.75
# yukawa_N_min            0.5
# yukawa_N_max            3.0
# yukawa_N_step           0.1

#if set , the yukawa coupling will be multiplied by sqrt(2 kappa), or lambda_N by (4 kappa^2) , or lambda_6_N by (8 kappa^3)
interpret_yukawa_as_continuum    0
interpret_lambda_as_continuum    0
interpret_lambda_6_as_continuum  0

scan_mag                    1 
mag                         0.5                        
mag_min                     0.01 
mag_max                     5
mag_step                    0.01

scan_stag                    1 
stag                         0.5                        
stag_min                     0.01 
stag_max                     5
stag_step                    0.01


#may use: [Ls]->Lxx [Lt]->Txx [k]->k_xxxx [l]->l_xxxx [l6]->l6_xxxx [y]->y_xxxx [m]->mag_xxxx [s]->stag_xxxx for the scanables: e.g. k_xx_yy with max and min)
#         [loop]->withLoop or noLoop    [bc]-pBC or aBC for periodic or antiperiodic bc   [det]->withDet or noNet, if use_improved_treeLevel is on or off
#         [1st] ->with1st or no1st if use_improved_firstOrder is on or off
#
outputfile                  output/potential/output_[Ls][Lt]_[y]_[l]_[l6]_[k]_[m]_[s].txt
