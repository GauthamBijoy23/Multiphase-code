#???# ====================================================
# Compiler setup
# ====================================================
FC = mpif90                # MPI Fortran compiler wrapper (gfortran backend)
#FC = /opt/nvidia/hpc_sdk/Linux_x86_64/24.5/comm_libs/mpi/bin/mpif90  (mod2)

# GPU architecture (70=V100, 80=A100, 90=H100)
GPU_ARCH = 70

#Usage: For switching betwwn 2 and 5 phase
# Use: make PHASE=5   → compiles 5-phase version flag -cpp
#      make PHASE=2   → compiles 2-phase version (default)
PHASE ?= 2

# Usage: For activating inviscid fluid
#   make AMUP=CONST   → amup is constant
#   make AMUP=VAR     → amup is variable (default)
AMUP ?= CONST
#AMUP will be activated only for the 5 phase model for casting problem only no need of K-Epsilon model with this

# Usage: For K-Epsilon model
#   make KE=OFF   → KE is not activated (default)
#   make KE=ON    → KE in ativated 
# Default KE turbulence = OFF
KE ?= OFF
# Right now available only for jameson

# Usage: For condensation problem 
#   make NP=OFF   → NP is not activated (default)
#   make NP=ON    → NP in ativated 
# Default NP = OFF
NP ?= OFF

# Usage: For transport of species 
#   make YP=OFF   → YP is not activated (default)
#   make YP=ON    → YP in ativated 
# Default YP = OFF
YP ?= OFF
# Right now available only for jameson

# Usage: For activation of chemical reactions:
#   make IREACT=OFF   → NP is not activated (default)
#   make IREACT=ON    → NP in ativated 
# Default IREACT = OFF
IREACT ?= OFF

# Usage: For activation of shock tube mode:
#   make ISHOCK=OFF   → Standard simulation (default input file)
#   make ISHOCK=ON    → Shock tube simulation (reads shock tube input file)
# Default ISHOCK = OFF
ISHOCKTUBE ?= OFF


#Flags for GPU compilation 
#FFLAGS = -O3 -acc  -Minfo=accel -mcmodel=medium -gpu=cc$(GPU_ARCH) -Wa,--noexecstack -cpp   (mod2)
FFLAGS = -O0 -g -fcheck=all -fbacktrace -mcmodel=medium -cpp

#Flags for CPU compilation
#FFLAGS = -O3 -mcmodel=medium  -Wa,--noexecstack -cpp

#For 2/5 phase
ifeq ($(PHASE),5)
    FFLAGS += -DFIVE_PHASE
endif

#For AMUP constant/variable
ifeq ($(AMUP),CONST)
    FFLAGS += -DAMUP_CONSTANT
endif

#For K-Epsilon turbulence model
ifeq ($(KE),ON)
    FFLAGS += -DKE_TURB
endif

#For condensation model
ifeq ($(NP),ON)
    FFLAGS += -DNP_P
endif

#For species transport model
ifeq ($(YP),ON)
    FFLAGS += -DYP_P
endif

#For chemical reaction 
ifeq ($(IREACT),ON)
    FFLAGS += -DIREACT
endif

#Shock tube tut
ifeq ($(ISHOCK),ON)
  FFLAGS += -DISHOCKTUBE
endif


# ====================================================
# Directories
# ====================================================
SRC_DIR = src
BUILD_DIR = build
EXEC = a.out

# List ALL your Fortran source files here:
SRC = \
  $(SRC_DIR)/mod_GlobalVariables.F90 \
  $(SRC_DIR)/mod_Fluxmodule.F90 \
  $(SRC_DIR)/mod_MatrixOps.F90 \
  $(SRC_DIR)/main.F90 \
  $(SRC_DIR)/sr_geometry.F90 \
  $(SRC_DIR)/sr_derivative.F90 \
  $(SRC_DIR)/sr_dissipationqrs.F90 \
  $(SRC_DIR)/sr_deallocate_variable.F90 \
  $(SRC_DIR)/sr_hllcbatten.F90 \
  $(SRC_DIR)/sr_hllczien.F90 \
  $(SRC_DIR)/sr_init_cond.F90 \
  $(SRC_DIR)/sr_jameson.F90 \
  $(SRC_DIR)/sr_jamesonqrs.F90 \
  $(SRC_DIR)/sr_mpiexchangecn.F90 \
  $(SRC_DIR)/sr_mpiexchangecu.F90 \
  $(SRC_DIR)/sr_mpiexchangedel.F90 \
  $(SRC_DIR)/sr_primitive.F90 \
  $(SRC_DIR)/sr_primitiveplusbc.F90 \
  $(SRC_DIR)/sr_Read_flow.F90 \
  $(SRC_DIR)/sr_RK.F90 \
  $(SRC_DIR)/sr_smooth.F90 \
  $(SRC_DIR)/sr_source.F90 \
  $(SRC_DIR)/sr_sourceqrs.F90\
  $(SRC_DIR)/sr_viscousgp.F90 \
  $(SRC_DIR)/sr_VK.F90

# Object files
OBJ = $(SRC:$(SRC_DIR)/%.F90=$(BUILD_DIR)/%.o)


# ====================================================
# Default target
# ====================================================
all: $(EXEC)

# ====================================================
# TIOGA target  (mod3)
# ====================================================
#tioga_build:
#	cd $(HOME)/Multiphase_code/HydraPhase_60_FPC_5mill/tioga/build && make


# ====================================================
# Link object files to executable
# ====================================================
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

# ====================================================
# Compile modules and other source files
# -J$(BUILD_DIR) tells gfortran where to place .mod files
# ====================================================
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.F90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -J$(BUILD_DIR) -o $@
        #$(FC) $(FFLAGS) -c $< -J$(BUILD_DIR) -o $@ (mod 2)


# ====================================================
# Create build directory if it doesn't exist
# ====================================================
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# ====================================================
# Clean object, module, and executable files
# ====================================================
clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/*.mod $(EXEC)

# ====================================================
# Phony targets
# ====================================================
.PHONY: all clean #tioga_build (mod3)

