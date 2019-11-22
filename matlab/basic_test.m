% basic test of randutv matlab interface. Barnett 11/13/19
clear
m = 10; n = 8;
A = randn(m,n);
[U,T,V] = randutv(A);
fprintf('||A - UTV^T|| = %.3g (should be close to emach)\n',norm(A - U*T*V'))
