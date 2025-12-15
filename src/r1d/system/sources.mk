r1d_system_=\
  r1d/system/System.cpp \
  r1d/system/SystemAccess.cpp \
  r1d/system/HomogeneousComparison.cpp 

r1d_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_system_:.cpp=.o))
