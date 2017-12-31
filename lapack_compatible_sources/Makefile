
# Defintions of variables.

CC      = gcc
CCFLAGS = -O -fopenmp
LD      = gcc
LDFLAGS = -O -fopenmp

# Defintions of rules.

simple_test.x : simple_test.o NoFLA_UTV_WY_blk_var2.o
	$(LD) $(LDFLAGS) \
            -o simple_test.x \
            simple_test.o \
            NoFLA_UTV_WY_blk_var2.o \
            /home/gquintan/mylib/lapack/liblapack_340_p4b64_gf.a \
            /home/gquintan/mylib/lapack/librefblas_340_p4b64_gf.a \
            -lm -lgfortran

            #### /usr/local/lapack/liblapack_340_p4b64_gf.a \
            #### /usr/local/lapack/mt_openblas/lib/libopenblas_haswellp-r0.2.14.a \

simple_test.o : simple_test.c
	$(CC) $(CCFLAGS) -O -fopenmp -c simple_test.c

NoFLA_UTV_WY_blk_var2.o : NoFLA_UTV_WY_blk_var2.c
	$(CC) $(CCFLAGS) -O -fopenmp -c NoFLA_UTV_WY_blk_var2.c

clean: 
	rm -f a.out *.x *.o *~ core

