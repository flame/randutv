# makefile for MEX interface to (lapack-compat version of) randutv
# Barnett 11/12/19

MWRAP = /usr/bin/mwrap
MEX = mex
# use MATLAB's shipped lin alg:
LIBS = -lmwblas -lmwlapack
#LIBS = -lblas -llapack
#LIBS = /usr/lib/openblas-base/liblapack.so /usr/lib/openblas-base/libblas.so

MEXFLAGS = -lgfortran -largeArrayDims
CC = gcc

# upon calling from matlab, get "Intel MKL ERROR"s, indicating matlab's blas etc
# being called, not openblas.
# See:
# https://www.mathworks.com/matlabcentral/answers/268503-linking-mex-with-lapack
# https://stackoverflow.com/questions/49794790/mkl-error-dgemm-wrong-entry
# https://www.mathworks.com/matlabcentral/answers/388042-using-external-lapack-versions-in-mex

# Try to preload blas+lapack:
# export LD_PRELOAD=/usr/lib/openblas-base/liblapack.so:/usr/lib/openblas-base/libblas.so:$LD_PRELOAD
#echo $LD_PRELOAD
#/usr/lib/openblas-base/liblapack.so:/usr/lib/openblas-base/libblas.so:/usr/lib/x86_64-linux-gnu/libgomp.so.1
# causes segfault crash in MATLAB.

SRC = ../lapack_compatible_sources/NoFLA_UTV_WY_blk_var2.c
# build its .o in this dir:
OBJ = NoFLA_UTV_WY_blk_var2.o

default: mexfile

# experts rebuild interface via mwrap
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
	$(MEX) $(CCFLAGS) gateway.c $(OBJ) $(MEXFLAGS) -lm -lgomp $(LIBS)

clean:
	rm -f *.mex* *.o

expertclean: clean
	rm -f randutv.m gateway.c
