#-----------------------------------------------------------------------
# Include *.mk makefile fragments

include ../../config.mk
include $(SRC_DIR)/prdc/sources.mk
include $(SRC_DIR)/prdc/patterns.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(prdc_cl_OBJS) 

clean:
	rm -f $(prdc_cl_OBJS) $(prdc_cl_OBJS:.o=.d) 

#-----------------------------------------------------------------------
# Include dependency files

-include $(prdc_cl_OBJS:.o=.d)
