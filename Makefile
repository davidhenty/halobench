MF=	Makefile

FC=	ftn
FFLAGS=

LFLAGS=

EXE=	halobench

SRC= \
	halobench.f90 \
	haloswap.f90 \
	benchclock.f90


#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

halobench.o: haloswap.o benchclock.o

clean:
	rm -f $(OBJ) $(EXE) core

tar:
	tar --exclude-vcs -cvf $(EXE).tar $(MF) $(SRC) benchio.pbs \
		defstriped/README striped/README unstriped/README
