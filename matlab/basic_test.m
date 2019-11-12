clear
m = 10; n = 8;
A = randn(m,n);
[U,T,V] = randutv(A);
norm(A - U*T*V)
