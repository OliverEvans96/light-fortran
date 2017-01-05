# Kelp light model Makefile

# Makefile links:
# http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
# http://mrbook.org/blog/tutorials/make/
# http://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules

# Project directories
BASE=/home/oliver/academic/research/kelp/fortran
BIN=$(BASE)/bin
SRC=$(BASE)/src
INC=$(BASE)/include

# Fortran Compiler
FC=gfortran

# Fortran Compilation flags
# Object files (.o)
OFLAGS=-J$(INC) -I$(INC) -c
# Binary files (executable)
BFLAGS=-J$(INC) -I$(INC)

all: test_interp

test_interp: rte_utils.o
	 $(FC) $(BFLAGS) $(SRC)/test_interp.f90 $(INC)/rte_utils.o -o $(BIN)/test_interp 

rte_utils.o:
	$(FC) $(OFLAGS) $(SRC)/rte_utils.f90 -o $(INC)/rte_utils.o

clean:
	rm -f $(INC)/*.mod $(INC)/*.o $(BIN)/*

ls:
	ls $(SRC) $(BIN) $(INC)
