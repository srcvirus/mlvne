CPLEX_ROOT = /opt/ibm/ILOG/CPLEX_Studio125
LIB_PATHS = -L$(CPLEX_ROOT)/cplex/lib/x86-64_sles10_4.1/static_pic -L$(CPLEX_ROOT)/concert/lib/x86-64_sles10_4.1/static_pic

INCLUDE_PATHS = -I$(CPLEX_ROOT)/cplex/include -I$(CPLEX_ROOT)/concert/include
LIBS = -lilocplex -lconcert -lcplex -lm -lpthread  -DIL_STD

FILES = mlvn.cc util.cc cplex_solver.cc vne_solution_builder.cc

all:
	g++ -std=c++0x -O3 $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o ml_vne

dbg:
	g++ -std=c++0x -DDBG -g $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o ml_vne
debug:
	g++ -std=c++0x -g $(LIB_PATHS) $(INCLUDE_PATHS) $(FILES) $(LIBS) -o ml_vne
