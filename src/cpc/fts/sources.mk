cpc_fts_= \
  cpc/fts/Simulator.cpp
  
cpc_fts_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_fts_:.cpp=.o))
