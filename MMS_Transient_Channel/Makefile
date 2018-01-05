include ../Makefile.site
include ../Makefile.common

all:
	(cd ForwardProblem; make forward_problem_TransientChannel.exe)
	(cd POD/SingleParameterSpace; make all)
	(cd InverseProblem; make all)
	(cd ODECO; make all)

clean:
	(cd ForwardProblem; make clean)
	(cd POD/SingleParameterSpace; make clean)
	(cd InverseProblem; make clean)
	(cd ODECO; make clean)	

spotless:
	(cd ForwardProblem; make spotless)
	(cd POD/SingleParameterSpace; make spotless)
	(cd InverseProblem; make spotless)
	(cd ODECO; make spotless)	
