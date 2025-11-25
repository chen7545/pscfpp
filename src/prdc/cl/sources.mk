prdc_cl_= \
  prdc/cl/Interaction.cpp

prdc_cl_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_cl_:.cpp=.o))

