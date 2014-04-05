function [] = test_omp_thm2()
% test_OMP_Thm2 - Driver routine to test signal recovery algorithm.
%                 Specifically, test that OMP works using the value of N
%                 given in Theorem 2 of Tropp 2007      
%
% Syntax: 
%  [] = test_omp_thm2()
%
% Inputs:
%  None
%
% Outputs:
%  None
%
% Example:
%  test_omp_thm2()
%
% Dependencies:
%  test_recovery
%  back_subs
% 
% TODO:
%  What is K from Thm 2 of Tropp 2007?
%
% Authors: JF,EY
% Revision history:
%  05 April 2014 - date written


%% For repeatability, set PRNG (Mersenne twister) and seed (seed = 0)
rng('default');


%% General parameters
d = 50; % signal length
delta = 0.1; % 0 < delta < 0.36, 1-2*delta <= OMP recovery probability
K = 5; % TODO what is a good value of K?


%% Generate reference signal and sparsify
m = ceil(0.95*d); % sparsity level
% reference signal
s_full = 2*rand([d 1])-1; % uniform distribution on [-1,1]
num_remove_inds = d-m;
perm_inds = randperm(d);
remove_inds = perm_inds(1:num_remove_inds);
sparse_inds = perm_inds(num_remove_inds+1:d);
s = zeros([d 1]);
s(sparse_inds) = s_full(sparse_inds); % sparse reference signal


%% Measurment vectors
N = ceil(K*m*log(d/delta)); % N from Thm 2 of Tropp 2007

mu_Phi = zeros([N d]); % mean
mu_Sigma = 1/N*eye([d d]); % covariance
Phi = mvnrnd(mu_Phi,mu_Sigma); % measurement matrix, columns are vectors from
                               % multivariate N(0,1/N)
v = Phi*s; % data vector


%% Perform OMP

% 1 - initialization
resid = v; % residual vector
Lambda = []; % index set (column vector)
Phi_t = []; % matrix of atoms (will be filled by OMP)

Q_Phi = [];
R_Phi = [];

for t = 1:m
   
   % 2 - find index
   % if there are multiple indexes, max returns the first index
   [maxval,lambda] = max(abs(transpose(resid)*Phi));

   % 3 - augment index set and matrix of atoms
   Lambda = cat(1, Lambda, lambda); % add index column
   Phi_t = cat(2, Phi_t, Phi(:,lambda)); % add column on right

   % 4 - solve least squares problem
   % update QR decomp for new QR
   [Q_Phi, R_Phi] = qrinsert(Q_Phi,R_Phi,t,Phi(:,lambda),'col');
   rhs = transpose(Q_Phi)*v;
   x = back_subs(R_Phi,rhs);
   
   %[q,r] = qr(Phi_t);
   %norm(Q_Phi-q,'fro')
   %norm(R_Phi-r,'fro')
   %norm(x-Phi_t\v)
  
   % 5 - update approximation and residual
   a = Phi_t*x;
   resid = v - a;
   
end

% 6 - construct approximate signal
s_hat = zeros([d 1]);
s_hat(Lambda) = x;

test_recovery(s,s_hat)

end % test_omp_thm2


%% Subroutines

function [x] = back_subs(R,b)
% back_subs - perform backward substitution on an upper trapezoidal matrix
%             this algorithm assumes that R is ''tall and skinny'', meaning
%             that num_rows >= num_cols.
%
% Syntax: 
%  [x] = back_subs(R,b)
%
% Inputs:
%  R - upper trapezoidal (triangular) matrix (e.g. from QR decomp)
%  b - RHS of linear system R*x=b
%
% Outputs:
%  x - solution of system
%
% Example:
%  A = rand([10 5]);
%  rhs = rand([10 1]);
%  [Q,R] = qr(A);
%  b = transpose(Q)*rhs;
%  x = back_subs(R,b);
%
% Dependencies:
%  None

ZERO_TOL = 10^(-14);

[m,n] = size(R);

if (m < n)
   error('back_subs: input matrix R is not ``tall and skinny``');
end

x = zeros([n 1]);

if (abs(R(n,n)) < ZERO_TOL)
   error('back_subs: input matrix R is singular to precision ZERO_TOL');
end

x(n) = b(n)/R(n,n);

for i=n-1:-1:1
   if (abs(R(i,i)) < ZERO_TOL)
      error('back_subs: input matrix R is singular to precision ZERO_TOL');
   end

   x(i) = (b(i) - R(i,i+1:n)*x(i+1:n))/R(i,i);
end

end % back_subs



function [] = test_recovery(s,s_hat)
% test_recovery - plot the reference signal, the approximation,
%                 and their absolute difference
%
% Syntax: 
%  [] = test_recovery(s,v)
%
% Inputs:
%  s - reference signal; d vector
%  s_hat - sparsified signal; d vector
%
% Outputs:
%  None
%
% Example:
%  test_recovery([1;2;3;4],[1.01;1.97;3.02;4.01]);
%
% Dependencies:
%  None


% plot ref. signal, approximation signal, and difference
subplot(3,1,1);
plot(s)
ylabel('ref sig')

subplot(3,1,2);
plot(s_hat)
ylabel('approx sig')

subplot(3,1,3);
plot(abs(s-s_hat))
ylabel('abs diff')
xlabel('index')

end % test_sparsify

