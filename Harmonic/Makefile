# Nome dell'eseguibile
EXEC = simulation

# Compilatore
FC = gfortran

# Flag di compilazione
# -I. serve per cercare i file .mod nella cartella corrente
# -g -fbacktrace sono utili per il debug se il codice crasha
FFLAGS = -O2 -Wall -I. -I/usr/local/include 

LDFLAGS = -lfftw3 -lm

# File sorgenti
SRC = kinds.f90 \
      fft_correlation_module.f90 \
      force_module.f90 \
      friction_module.f90 \
      main.f90

# File oggetto (sostituisce .f90 con .o)
OBJ = $(SRC:.f90=.o)

# Regola principale
all: $(EXEC)

# Link finale (Creazione eseguibile)
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compilazione generica dei file oggetto
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Definiscono l'ordine di compilazione corretto per i moduli Fortran.

kinds.o: kinds.f90

fft_correlation_module.o: fft_correlation_module.f90 kinds.o

force_module.o: force_module.f90 kinds.o

friction_module.o: friction_module.f90 fft_correlation_module.o kinds.o

main.o: main.f90 force_module.o friction_module.o kinds.o

# Pulizia
clean:
	rm -f *.o *.mod $(EXEC)

.PHONY: all clean
