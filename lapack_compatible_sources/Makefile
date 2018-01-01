#
# Definitions of variables for compiling.
#
LAPACK_LIB = /usr/local/mylib/lapack/liblapack_340_p4b64_gf.a
BLAS_LIB   = /usr/local/mylib/lapack/librefblas_340_p4b64_gf.a

CC         = gcc
CCFLAGS    = -O -fopenmp
LD         = gcc
LDFLAGS    = -O -fopenmp

#
# Definitions of files and rules.
#
OBJS = basic_test.o NoFLA_UTV_WY_blk_var2.o

basic_test.x : $(OBJS)
	$(LD) $(LDFLAGS) \
            -o $@ \
            $(OBJS) \
            $(LAPACK_LIB) \
            $(BLAS_LIB) \
            -lm -lgfortran

%.o : %.c
	$(CC) $(CCFLAGS) -c $< -o $@

clean: 
	rm -f a.out *.x *.o *~ core

