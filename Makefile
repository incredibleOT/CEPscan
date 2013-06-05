CXX=g++

# CXXFLAGS=-DF_
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic 

# LIBS=""

# INCL=""

OFILES = constrainedEffectivePotential.o CEPscan_helper.o


CEPFiles= constrainedEffectivePotential_computeFunctionOnly.cc constrainedEffectivePotential_computeGradientOnly.cc constrainedEffectivePotential_computeFunctionAndGradient.cc constrainedEffectivePotential_computeBosonicLoop.cc constrainedEffectivePotential_computeSecondDerivative.cc

INCL=-I/opt/products/gsl/1.15/include

LIBSLOC=-L/opt/products/gsl/1.15/lib64

LIBS= -lgsl -lgslcblas -lm -static
# ##########################################

CEPscan: CEPscan.o ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}

CEPscan.o: CEPscan.cc ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<

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

CEPscan_helper.o: CEPscan_helper.cc CEPscan_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 

clean:
	rm *.o *~

clean_exec:
	rm testCEPscan timingTestCEP

