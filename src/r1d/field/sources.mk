r1d_field_=\
  r1d/field/GeometryMode.cpp \
  r1d/field/Domain.cpp \
  r1d/field/FieldIo.cpp 

r1d_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_field_:.cpp=.o))

