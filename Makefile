# Kelp light model Makefile

##################
# Makefile links #
##################

# http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html
# http://mrbook.org/blog/tutorials/make/
# http://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules

#########
# Flags #
#########

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

###############
# Executables #
###############

all: test_interp test_rte2d test_vsf

test_interp: utils.o
	 $(FC) $(BFLAGS) $(SRC)/test_interp.f90 $(INC)/utils.o -o $(BIN)/test_interp 

test_rte2d: rte2d.o

test_vsf: rte_core.o
	 $(FC) $(BFLAGS) -g -fcheck=all -Wall $(SRC)/test_vsf.f90 $(INC)/rte_core.o $(INC)/utils.o -o $(BIN)/test_vsf 

################
# Object files #
################

rte2d.o: rte_core.o
	$(FC) $(OFLAGS) $(SRC)/rte2d.f90 -o $(INC)/rte2d.o

rte_core.o: utils.o
	$(FC) $(OFLAGS) $(SRC)/rte_core.f90 -o $(INC)/rte_core.o

utils.o:
	$(FC) $(OFLAGS) $(SRC)/utils.f90 -o $(INC)/utils.o

clean:
	rm -f $(INC)/*.mod $(INC)/*.o $(BIN)/*

ls:
	ls $(SRC) $(BIN) $(INC)
