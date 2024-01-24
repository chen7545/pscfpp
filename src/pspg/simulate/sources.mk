# Include source files lists from subdirectories
include $(SRC_DIR)/pspg/simulate/mcmove/sources.mk
include $(SRC_DIR)/pspg/simulate/bdstep/sources.mk
include $(SRC_DIR)/pspg/simulate/trajectory/sources.mk
include $(SRC_DIR)/pspg/simulate/analyzer/sources.mk

pspg_simulate_= \
  $(pspg_simulate_mcmove_) \
  $(pspg_simulate_bdstep_) \
  $(pspg_trajectory_) \
  $(pspg_simulate_analyzer_) \
  pspg/simulate/Simulator.cu \
  pspg/simulate/SimulatorFactory.cu
  
pspg_simulate_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_simulate_:.cu=.o))

