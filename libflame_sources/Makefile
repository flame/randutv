
# Defintions of variables.

LIBFLAME_PATH = /usr/local/mylib/libflame-pc34-r11104-gc

CC      = gcc
CCFLAGS = -O -fopenmp -I$(LIBFLAME_PATH)/include-x86_64-r -msse3
LD      = gcc
LDFLAGS = -O -fopenmp
LDFLAGS = -O -fopenmp

# Defintions of rules.

simple_test.x : simple_test.o FLA_UTV_UT_blk_var1.o
	$(LD) $(LDFLAGS) \
            -o simple_test.x \
            simple_test.o \
            FLA_UTV_UT_blk_var1.o \
            -L$(LIBFLAME_PATH)/lib -lflame-x86_64-r \
            /usr/local/mylib/lapack/liblapack_340_p4b64_gf.a \
            /usr/local/mylib/lapack/librefblas_340_p4b64_gf.a \
            -lgfortran \
            -lm

%.o : %.c
	$(CC) $(CCFLAGS) -c $< -o $@

clean: 
	rm -f a.out *.x *.o *~ core

