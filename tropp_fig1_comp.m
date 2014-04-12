function [] = tropp_fig1_comp(matfile)
% tropp_fig1_comp - routine to reproduce data for figure 1 on Tropp 2007
%
% Syntax: 
%  [] = tropp_fig1_comp()
%  [] = tropp_fig1_comp(matfile)
%
% Inputs:
%  matfile - save data to matfile (it defaults to a different value for 
%              linux and mac
%
% Examples:
%  >> tropp_fig1_comp() % use default save file
%  >> tropp_fig1_comp('savefile.mat') % use the given save file
%
% Outputs:
%  None
%
% Dependencies:
%  omp_alg
%  check_recovery
% 
% TODO:
%  Check if default save file works for mac
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

%% default savefile
if (nargin == 0)
   if strcmp(computer(),'GLNXA64')
      matfile = 'fig1_data_glnx64.mat';
   elseif strcmp(computer(),'MACI64')
      matfile = 'fig1_data_mac64.mat';
   else
      error('tropp_fig1_comp: default file not given for the arch used');
   end
end

save(matfile);


%% figure parameters
num_sigs = 50; % TODO change this to 1000
m_vec = [4,12,20,28,36];
N_vec = 1:20:d; % TODO finer resolution
percent_recovered = zeros([numel(N_vec) numel(m_vec)]);


%% generate signals and try to recovery them
for i_m = 1:numel(m_vec)
   fprintf(1,'m = %d  N = ',m_vec(i_m));
   for i_N = 1:numel(N_vec)
      fprintf(1,'%d ',N_vec(i_N));
      mu_Phi = zeros([N_vec(i_N) d]); % mean
      Sigma_Phi = 1/N_vec(i_N)*eye([d d]); % covariance
      Phi = mvnrnd(mu_Phi,Sigma_Phi); % measurement matrix, columns

      for sig_ind = 1:num_sigs
         s = gen_sig(d,m_vec(i_m));
         v = Phi*s;
         s_hat = omp_alg(m_vec(i_m),Phi,v);
         
         percent_recovered(i_N, i_m) = percent_recovered(i_N, i_m) + ...
            check_recovery(s, s_hat);
      end
   end
   fprintf(1,'\n');
end

percent_recovered = percent_recovered / num_sigs * 100;

save(matfile);

end % tropp_fig1_comp
