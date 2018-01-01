#
# Definitions of variables for compiling.
#
LIBFLAME_PATH = /usr/local/mylib/libflame-pc34-r11104-gc
LAPACK_PATH   = /usr/local/mylib/lapack
BLAS_PATH     = /usr/local/mylib/lapack

CC      = gcc
CCFLAGS = -O -fopenmp -I$(LIBFLAME_PATH)/include-x86_64-r
LD      = gcc
LDFLAGS = -O -fopenmp
LDFLAGS = -O -fopenmp

#
# Definitions of files and rules.
#
OBJS = basic_test.o FLA_UTV_UT_blk_var1.o

basic_test.x : $(OBJS)
	$(LD) $(LDFLAGS) \
            -o $@ \
            $(OBJS) \
            -L$(LIBFLAME_PATH)/lib -lflame-x86_64-r \
            $(LAPACK_PATH)/liblapack_340_p4b64_gf.a \
            $(BLAS_PATH)/librefblas_340_p4b64_gf.a \
            -lgfortran \
            -lm

%.o : %.c
	$(CC) $(CCFLAGS) -c $< -o $@

clean: 
	rm -f a.out *.x *.o *~ core

