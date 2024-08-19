# Define the Fortran compiler
FC = gfortran

# Compiler flags
FFLAGS = -O2 -Wall

# The target executable
TARGET = SpinCross.x

# The directory containing the source files
SRC_DIR = src

# List of Fortran source files
SRC_FILES = $(SRC_DIR)/EffectieGrad.f90 \
            $(SRC_DIR)/General_QuantChem_Interface.f90 \
            $(SRC_DIR)/InitOptimizer.f90 \
            $(SRC_DIR)/InputReader.f90 \
            $(SRC_DIR)/main_driver.f90 \
            $(SRC_DIR)/mainvar.f90 \
            $(SRC_DIR)/Optimizers.f90 \
            $(SRC_DIR)/Orca_Interface.f90 \
            $(SRC_DIR)/Utilities.f90

# Object files
OBJ_FILES = $(SRC_FILES:.f90=.o)

# Default target: compile the program
all: $(TARGET)

# Rule to build the target executable
$(TARGET): $(OBJ_FILES)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ_FILES)

# Rule to compile each Fortran source file into an object file
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up object files and the executable
clean:
	rm -f $(OBJ_FILES) $(TARGET)

# Phony targets
.PHONY: all clean

