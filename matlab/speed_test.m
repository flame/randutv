% speed test of randutv matlab interface vs full svd.
% real matrices only (exists complex version?)
% Alex Barnett 11/13/19
clear
%maxNumCompThreads(12);   % choose whatever you like
n = 5e3;
A = rand(n,n)-0.5;
disp('starting randutv...')
tic;
[U,T,V] = randutv(A);           % my default
%[U,T,V] = randutv(A,96,4,1);   % more control (64 beats 32 or 128)
t1=toc; fprintf('randutv n=%d: %.3g s\n',n,t1)
%fprintf('||A - UTV^T|| = %.3g (should be close to emach)\n',norm(A - U*T*V'))
tic; [UU,S,VV] = svd(A);
t2=toc; fprintf('full svd: %.3g s\n',t2)
%fprintf('||A - USV^T|| = %.3g (should be close to emach)\n',norm(A - UU*S*VV'))
fprintf('speedup factor %.3g\n',t2/t1)
fprintf('max ratio of diag el of T to sing val: %.3g\n',max(diag(T)./diag(S)))

% I see for n=1e4: 2x speedup (& max ratio < 1.8) on matlab R2019a, xeon 12-core
