pscf_cpu_= \
  pscf/cpu/complex.cpp \
  pscf/cpu/VecOpCx.cpp 

pscf_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cpu_:.cpp=.o))

