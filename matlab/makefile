# makefile for MEX interface to (lapack-compat version of) randutv
# Barnett 11/12/19

MEX = mex
# use MATLAB's shipped lin alg (MKL, needing new ptr sizes below)
LIBS = -lmwblas -lmwlapack
# (there is no way to avoid linking to these, in my experience)

MEXFLAGS = -v -largeArrayDims
# for recent matlab versions >= R2017b
#MEXFLAGS = -v -R2017b

CC = gcc

# https://www.mathworks.com/help/matlab/matlab_external/calling-lapack-and-blas-functions-from-mex-files.html
# describes that need to use ptrdiff_t type for int ptrs going to Matlab's MKL
# blas/lapack. This entails hacking the .c and .h to use long int.

# hacked local copy using long int everywhere instead of int...
SRC = NoFLA_UTV_WY_blk_var2_int8.c
OBJ = NoFLA_UTV_WY_blk_var2_int8.o

default: mexfile

# experts rebuild interface via mwrap
MWRAP = mwrap
expert: expertclean gateway.c mexfile

$(OBJ): $(SRC)
	$(CC) $(SRC) -c -o $(OBJ) -O -fPIC

# Tell mwrap to make the interface and the name of the final MEX file:
# (note -mb here makes separate .m file for each @ function in the .mw)
gateway.c: randutv.mw makefile
	$(MWRAP) -list -mex gateway -mb randutv.mw
	$(MWRAP) -mex gateway -c gateway.c randutv.mw

# How to make MATLAB/MEX compile the interface (whatever libs needed):
mexfile: $(OBJ) makefile
	$(MEX) $(MEXFLAGS) gateway.c $(OBJ) -lm -lgomp $(LIBS)

clean:
	rm -f *.mex* *.o

expertclean: clean
	rm -f randutv.m gateway.c
