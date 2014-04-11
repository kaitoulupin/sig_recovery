function [] = tropp_fig1()
% tropp_fig1 - routine to reproduce figure 1 on Tropp 2007
%
% Syntax: 
%  [] = tropp_fig1()
%
% Inputs:
%  None
%
% Outputs:
%  None
%
% Dependencies:
%  omp_alg
%  check_recovery
% 
% TODO:
%
%
% Authors: JF,EY
% Revision history:
%  11 April 2014 - date written
                    

%% For repeatability, set PRNG (Mersenne twister) and seed (seed = 0)
rng('default');


%% General parameters
d = 256; % signal length
delta = 0.1; % 0 < delta < 0.36, 1-2*delta <= OMP recovery probability
K = 5; %For our particular choices, K<=6.7874 is good.  See GC.pdf


%% figure parameters
num_sigs = 100; % TODO change this to 1000
m_vec = [4;12;20;28;36];
N_vec = 1:10:d;
percent_recoverd = zeros([numel(N_vec) numel(m_vec)]);


%% generate signals and try to recovery them
for m = m_vec
   for N = N_vec
      mu_Phi = zeros([N d]); % mean
      Sigma_Phi = 1/N*eye([d d]); % covariance
      Phi = mvnrnd(mu_Phi,Sigma_Phi); % measurement matrix, columns

      for sig_ind = 1:num_sigs
         s = gen_sig(d,m);
         error('not implemented');
      end
   end
end



%%% Generate reference signal and sparsify
%m = ceil((1-0.95)*d); % sparsity level
%s = gen_sig(d,m);
%
%%% Measurment vectors
%%N = ceil(K*m*log(d/delta)); % N from Thm 2 of Tropp 2007
%N = 150;
%
%mu_Phi = zeros([N d]); % mean
%mu_Sigma = 1/N*eye([d d]); % covariance
%Phi = mvnrnd(mu_Phi,mu_Sigma); % measurement matrix, columns are vectors from
%                               % multivariate N(0,1/N)
%v = Phi*s; % data vector
%size(Phi)
%
%%% Perform OMP
%[s_hat] = omp_alg(m,Phi,v);
%
%% 7 - check to see if s is recovered
%check_recovery(s,s_hat)



end % test_omp_smallN

