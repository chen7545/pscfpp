pspc_mcmove_= \
  pspc/mcmove/McMove.cpp \
  pspc/mcmove/McMoveFactory.cpp \
  pspc/mcmove/McState.cpp \
  pspc/mcmove/McSimulator.cpp 

pspc_mcmove_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_mcmove_))
pspc_mcmove_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_mcmove_:.cpp=.o))

