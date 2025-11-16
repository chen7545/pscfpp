cpc_field_= \
  cpc/field/FieldIo.cpp \
  cpc/field/Domain.cpp 

cpc_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_field_:.cpp=.o))

