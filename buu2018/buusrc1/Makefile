#-------------- Adjust the following variables as necessary -------------
#
# name of executable:
exec = BUU.x
# path of *.o files (~ for home directory doesn't work!?):
objpath = ../buuobj/

# name of FORTRAN compiler:
FC = gfortran
#FC = gfortran -pipe
#FC = xlf
# flags for FORTRAN compiler:
#FFLAGS = -ffpe-trap=zero -Wall -fbounds-check -fd-lines-as-comments -mcmodel medium -shared-intel
FFLAGS = -ffpe-trap=zero -Wall -fbounds-check -fd-lines-as-comments
#FFLAGS = -O3 -qstrict -qarch=auto -qtune=auto -qcache=auto
#FFLAGS = -ffpe-trap=zero,overflow,underflow -Wall -fbounds-check -fd-lines-as-comments


# name of c compiler:
CC = cc
#
#------------------------------------------------------------------------


OBJ = $(patsubst %.f,$(objpath)%.o,$(wildcard *.f)) $(objpath)cptime.o $(objpath)f77flush.o

$(exec): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $(exec)

$(objpath)cptime.o: cptime.c
	-mkdir -p $(objpath)
	$(CC) -c $(CFLAGS) cptime.c -o $(objpath)cptime.o

$(objpath)f77flush.o: f77flush.c
	-mkdir -p $(objpath)
	$(CC) -c $(CFLAGS) f77flush.c -o $(objpath)f77flush.o

$(objpath)%.o: %.f common cominput coucom com_kminu commsp commonthreef com_cont_epair com_pert
	-mkdir -p $(objpath)
	$(FC) -c $(FFLAGS) $< -o $@

.PHONY: clean

clean:
	-rm $(OBJ) *~ \#*\#
	-rmdir $(objpath)
