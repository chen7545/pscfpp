#-----------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/rpc/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(rpc_scft_sweep_OBJS) 

clean:
	rm -f $(rpc_scft_sweep_OBJS) $(rpc_scft_sweep_OBJS:.o=.d)

veryclean:
	$(MAKE) clean
	-rm -f *.o 
	-rm -f *.d 

#-----------------------------------------------------------------------
# Include dependency files

-include $(rpc_OBJS:.o=.d)
-include $(prdc_OBJS:.o=.d)
-include $(pscf_OBJS:.o=.d)
-include $(util_OBJS:.o=.d)
