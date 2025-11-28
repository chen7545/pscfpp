rpg_field_= \
  rpg/field/HostDArrayComplex.cu \
  rpg/field/FieldIo.cu \
  rpg/field/Domain.cu \
  rpg/field/WFields.cu \
  rpg/field/CFields.cu \
  rpg/field/Mask.cu

rpg_field_OBJS=\
	  $(addprefix $(BLD_DIR)/, $(rpg_field_:.cu=.o))

