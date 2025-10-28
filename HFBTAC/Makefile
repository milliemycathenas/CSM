mkllib   := /opt/intel/mkl/lib/intel64
mklinc   := /opt/intel/mkl/include

TAR      := ~/programs/HFBTAC/HFBTAC
BUILD    := build
MODULE_SRC := modules.f90
OTHER_SRCS := $(filter-out $(MODULE_SRC),$(wildcard *.f90))
SRCS       := $(MODULE_SRC) $(OTHER_SRCS)
MODULE_OBJ := $(BUILD)/module.o
MODULE_MOD := $(BUILD)/module.mod
OTHER_OBJS := $(patsubst %.f90,$(BUILD)/%.o,$(OTHER_SRCS))
OBJS       := $(MODULE_OBJ) $(OTHER_OBJS)

FCC     := ifx
LD_FLAGS := -I$(mklinc) -L$(mkllib) \
            -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
            -lmkl_lapack95_lp64 -liomp5 -lpthread

.PHONY: all clean

all: $(BUILD) $(TAR)

$(TAR): $(OBJS)
	$(FCC) -o $@ $^ $(LD_FLAGS) -g -O0

$(MODULE_MOD) $(MODULE_OBJ): $(MODULE_SRC) | $(BUILD)
	$(FCC) -c $< \
	    -module $(BUILD) -I$(BUILD) -I$(mklinc) \
	    -g -O0 -o $(MODULE_OBJ)

$(BUILD)/%.o: %.f90 | $(BUILD) $(MODULE_MOD)
	$(FCC) -c $< \
	    -module $(BUILD) -I$(BUILD) -I$(mklinc) \
	    -g -O0 -o $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD) $(TAR)
