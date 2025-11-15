cpc_solvers_= \
  cpc/solvers/Propagator.cpp \
  cpc/solvers/Block.cpp

cpc_solvers_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_solvers_:.cpp=.o))

