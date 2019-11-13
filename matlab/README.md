MATLAB interface to randutv, via the lapack-compatible version.

(Uses mwrap by David Bindel, although user does not need to know this.)

Alex Barnett 11/13/19

From this directory, check in the `makefile` the choice of compiler, etc,
then
```
make
```
to compile and build the MEX file.
From in MATLAB, run
```
basic_test
```
to check it worked. Run
```
speed_test
```
to compare randUTV against MATLAB's full SVD.

Here is the documentation:
```
 RANDUTV  compute UTV decomposition (approximate SVD) via randomized methods

 [U,T,V] = randutv(A) computes a rank-revealing UTV decomposition such that
   A = U*T*V', with U and V orthogonal, and T upper-triangular. This is
   designed to be a faster replacement for the SVD, with nearly the same
   accuracy.  A must be a double-precision real matrix.

 [U,T,V] = randutv(A,nb,pp,niter) also controls options:
   nb:       Block size. Usual values for nb are 32, 64, etc.
   pp:       Oversampling size. Usual values for pp are 5, 10, etc.
   niter:    Number of "power" iterations. Usual values are 0,1,2.

 Notes:
 1) Is a MEX interface to multi-threaded lapack-compatible implementation
    NoFLA_UTV_WY_blk_var2.c by Martinsson, Quintana-Orti and Heavner, 2016,
    which makes use of BLAS/LAPACK (those in matlab's MKL are used).
    For details see: https://github.com/flame/randutv
 2) A is copied to new matrix which is overwritten by T. Matlab interface
    could be extended so that input A is overwritten instead by T.
```

Tasks remaining: no-U,V version; complex version.






