prdc_rl_= \
  prdc/rl/fieldIoUtil.cpp

prdc_rl_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_rl_:.cpp=.o))

