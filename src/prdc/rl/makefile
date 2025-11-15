#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(prdc_rl_OBJS) 

clean:
	rm -f $(prdc_rl_OBJS) $(prdc_rl_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_rl_OBJS:.o=.d)
