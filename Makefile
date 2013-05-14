CXX=g++

# CXXFLAGS=-DF_
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic 

# LIBS=""

# INCL=""

OFILES = constrainedEffectivePotential.o


CEPFiles= constrainedEffectivePotential_computeFunctionOnly.cc constrainedEffectivePotential_computeGradientOnly.cc constrainedEffectivePotential_computeFunctionAndGradient.cc

INCL=-I/opt/products/gsl/1.15/include

LIBSLOC=-L/opt/products/gsl/1.15/lib64

LIBS= -lgsl -lgslcblas -lm -static
# ##########################################

testCEPscan: testCEPscan.o ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}

testCEPscan.o: testCEPscan.cc ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<

timingTestCEP: timingTestCEP.o ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}

timingTestCEP.o: timingTestCEP.cc ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<
###########################################3
###   OFILE   ###

constrainedEffectivePotential.o: constrainedEffectivePotential.cc constrainedEffectivePotential.h ${CEPFiles}
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 

clean:
	rm testCEPscan timingTestCEP *.o 

