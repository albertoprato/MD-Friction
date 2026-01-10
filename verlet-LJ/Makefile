# Executable name
EXEC = simulation

# Compiler
FC = gfortran

# Compilation flags
# -I. is needed to search for .mod files in the current directory
# -g -fbacktrace are useful for debugging if the code crashes
FFLAGS = -O2 -Wall -I. -I/usr/local/include

# Linker flags (FFTW3 library and Math library)
LDFLAGS = -lfftw3 -lm

# Source files
SRC = kinds.f90 \
      force_module.f90 \
      fft_correlation_module.f90 \
      minimization_module.f90 \
      velocity_init_module.f90 \
      friction_module.f90 \
      main.f90

# Object files (replaces .f90 with .o)
OBJ = $(SRC:.f90=.o)

# Main rule
all: $(EXEC)

# Final Link (Create executable)
$(EXEC): $(OBJ)
  $(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Generic compilation of object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# ==========================================
# Module Dependencies definition
# ==========================================

# Kinds module
kinds.o: kinds.f90

# Force module (depends on kinds)
force_module.o: force_module.f90 kinds.o

# FFT module (depends on kinds)
fft_correlation_module.o: fft_correlation_module.f90 kinds.o

# Velocity inizialitation module (depends on kinds)
velocity_init_module.o: kinds.o

# Optimization module (depends on kinds and force)
minimization_module.o: minimization_module.f90 force_module.o kinds.o

# Friction module (depends on kinds and FFT)
friction_module.o: friction_module.f90 fft_correlation_module.o kinds.o

# Main program (depends on ALL modules)
main.o: main.f90 force_module.o velocity_init_module.o minimization_module.o friction_module.o kinds.o

# Cleaning rule
clean:
	rm -f *.o *.mod $(EXEC)

.PHONY: all clean
