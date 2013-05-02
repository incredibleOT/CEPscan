CXX=g++

# CXXFLAGS=-DF_
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic 

# LIBS=""

# INCL=""

OFILES = constrainedEffectivePotential.o

effPotFiles=  effectivePotential_computePropagatorSum.cc effectivePotential_getFermionicContribution.cc effectivePotential.h

INCL=-I/opt/products/gsl/1.15/include
# ##########################################

testCEPscan: testCEPscan.o ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} -o $@ $^ ${LIBS}

testCEPscan.o: testCEPscan.cc ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<




###########################################3
###   OFILE   ###

constrainedEffectivePotential.o: constrainedEffectivePotential.cc constrainedEffectivePotential.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<
