pscf_cpu_= \
  pscf/cpu/complex.cpp \
  pscf/cpu/VecOpCx.cpp \
  pscf/cpu/ReduceCx.cpp 

pscf_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cpu_:.cpp=.o))

