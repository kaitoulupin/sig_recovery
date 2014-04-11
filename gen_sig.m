function [s] = gen_sig(d,m)
% gen_sig - generate a sparse random signal
%
% Syntax: 
%  [s] = check_recovery(d,m)
%
% Inputs:
%  d - length of signal
%  m - sparsity level (max number of nonzeros)
%
% Outputs:
%  s - signal
%
% Dependencies:
%  None
% 
% TODO:
%
% Authors: JF,EY
% Revision history:
%  11 April 2014 - date written

%% Generate reference signal and sparsify
s_full = 2*rand([d 1])-1; % uniform distribution on [-1,1]
num_remove_inds = d-m;
perm_inds = randperm(d);
remove_inds = perm_inds(1:num_remove_inds);
sparse_inds = perm_inds(num_remove_inds+1:d);
s = zeros([d 1]);
s(sparse_inds) = s_full(sparse_inds); % sparse reference signal

end % gen_sig
