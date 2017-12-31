
# Defintions of variables.

CC      = gcc
CCFLAGS = -O -fopenmp
LD      = gcc
LDFLAGS = -O -fopenmp

# Defintions of rules.

simple_test.x : simple_test.o
	$(LD) $(LDFLAGS) \
            -o simple_test.x \
            simple_test.o \
            NoFLA_UTV_WY_blk_var2.c \
            /usr/local/mylib/lapack/liblapack_340_p4b64_gf.a \
            /usr/local/mylib/lapack/librefblas_340_p4b64_gf.a \
            -lm -lgfortran

%.o : %.c
	$(CC) $(CCFLAGS) -c $< -o $@

clean: 
	rm -f a.out *.x *.o *~ core

