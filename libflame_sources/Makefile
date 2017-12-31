
# Defintions of variables.

LIBFLAME_PATH = /home/gquintan/mylib/libflame-pc34-r11104-gc

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
            /home/gquintan/mylib/lapack/liblapack_340_p4b64_gf.a \
            /home/gquintan/mylib/lapack/librefblas_340_p4b64_gf.a \
            -lgfortran \
            -lm

simple_test.o : simple_test.c
	$(CC) $(CCFLAGS) -c simple_test.c

FLA_UTV_UT_blk_var1.o : FLA_UTV_UT_blk_var1.c
	$(CC) $(CCFLAGS) -c FLA_UTV_UT_blk_var1.c

clean: 
	rm -f a.out *.x *.o *~ core

