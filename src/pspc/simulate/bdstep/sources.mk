pspc_simulate_bdstep_= \
  pspc/simulate/bdstep/BdSimulator.cpp \
  pspc/simulate/bdstep/BdStep.cpp \
  pspc/simulate/bdstep/BdStepFactory.cpp \
  pspc/simulate/bdstep/ExplicitBdStep.cpp \
  pspc/simulate/bdstep/MidstepBdStep.cpp \
  pspc/simulate/bdstep/PredCorrBdStep.cpp \
  pspc/simulate/bdstep/LMBdStep.cpp 
  
pspc_simulate_bdstep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_simulate_bdstep_:.cpp=.o))

