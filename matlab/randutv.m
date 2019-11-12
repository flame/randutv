function [U,T,V]=randutv(A,nb,pp,niter)
% RANDUTV  compute UTV decomposition (approximate SVD) via randomized methods
%
% [U,T,V] = randutv(A) basic usage
%   A must be a double-precision real matrix.
%
% [U,T,V] = randutv(A,nb,pp,niter) also controls options:
%   nb:       Block size. Usual values for nb are 32, 64, etc.
%   pp:       Oversampling size. Usual values for pp are 5, 10, etc.
%   niter:    Number of "power" iterations. Usual values are 0,1,2.
%
% Notes:
% 1) MEX interface to multi-threaded lapack-compatible implementation
%    NoFLA_UTV_WY_blk_var2.c by Martinsson, Quintana-Orti and Heavner, 2016.
%    For details see: https://github.com/flame/randutv
% 2) A is copied to new matrix which is overwritten by T. Matlab interface
%    could be extended so that input A is overwritten instead by T.

% Alex Barnett 11/12/19

if nargin<2, nb = 64; end    % defaults
if nargin<3, pp = 5; end
if nargin<4, niter = 1; end

if ~isa(A,'double') || ~isreal(A)
  error('A must be real double-prec array!');
end
[m,n] = size(A);

build_u = 1;    % basic task for now, output U,V
build_v = 1;

% documentation was unclear: A is overwritten by T! Thus need copy A to T:
%T = zeros(m,n);  % allocates new T
T = A;           % forces to copy A in

% MATLAB preallocation of U,V needed since it's 'inout' type, to allow
% later case of U,V being inputs too (which I don't understand...)
U = zeros(m,m);
V = zeros(n,n);

sizT = m*n; sizU = m*m; sizV = n*n;
% note all leading dims inputs are the # rows, ie, A,U,V, all contiguous.
mex_id_ = 'NoFLA_UTV_WY_blk_var2(i int, i int, io double[x], i int, i int, i int, i int, io double[x], i int, i int, i int, i int, io double[x], i int, i int, i int, i int)';
[T, U, V] = gateway(mex_id_, m, n, T, m, build_u, m, m, U, m, build_v, n, n, V, n, nb, pp, niter, sizT, sizU, sizV);
% here T,U,V were the outputs

U = reshape(U,[m,m]);
V = reshape(V,[n,n]);

