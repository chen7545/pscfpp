#-----------------------------------------------------------------------
# Source and object file lists for src/cpc 

# Include source list files from subdirectories
include $(SRC_DIR)/cpc/solvers/sources.mk

# List of source files in src/cpc
cpc_= \
  $(cpc_solvers_) 

# List of object file targets
cpc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_:.cpp=.o))

#-----------------------------------------------------------------------
# Path and rule for the cpc/libcpc.a library 

cpc_LIB=$(BLD_DIR)/cpc/libcpc.a

$(cpc_LIB): $(cpc_OBJS)
	$(AR) rcs $(cpc_LIB) $(cpc_OBJS)

