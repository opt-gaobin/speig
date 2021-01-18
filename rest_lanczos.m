function [ei,U,parout] = rest_lanczos(S,k,parin)
% This function computes some smallest eigenvalues of positive definite
% Hamiltonian matrices
%              H = JS
% where J is the Poisson matrix, S is an spd matrix. The approach is to
% used symplectic restarted Lanczos procedure.
% INPUT
%   S: the spd matrix such that H = JS
%   k: the number of smallest eigenvalues of interest
%   p: the addition Lanczos steps performed to improve the accuracy
%   vtil: the initial guess
%   tolcoeff: tolerance in compute the coefficients (see coef.m)
%   tolrest: tolerance for stopping restart (see lanczosing.m)
%   maxrest: maximal number of times to restart (see lanczosing.m)
% OUTPUT
%   ei: the (positive) imaginary part of the k smallest eigenvalues of H = JS
%   nrest: number of restarted performed
%   sdiff: to check the convergence (see lanczosing.m)
%   U*diag(Q1,Q1): symplectic Stiefel matrix in [1, Theorem 2]
%   where Q1 is the orthogonal matrix from running eig(T)
%   mess: information about the success
%   NOTE: 
% Reference
% [1] P. Amodio, On the computation of few eigenvalues of positive definite
% Hamiltonian matrices, Future Generation Computer Systems 22(2006) 403-411
% Author(s)
% NT Son, UCLouvain, 2020-02-12
% NT Son, 2020-03-28: computing the smallest 
% NT Son, 2020-12-17: computing the corresponding symplectic eigenvector
%% default setting
if nargin < 2; error('There must be at least two input'); end
if nargin < 3; parin = []; end
[Schol,chS] = chol(S,'lower');
if chS ~= 0
    error('S is not spd')
end
if ~isfield(parin, 'p');        parin.p = k; end
if ~isfield(parin, 'vtil');     parin.vtil = randn(size(S,1),1); end
if ~isfield(parin, 'tolc');     parin.tolc = 1e-6; end
if ~isfield(parin, 'tolr');     parin.tolr = 1e-7; end
if ~isfield(parin, 'maxr');     parin.maxr = 2*k; end
% cop
p = parin.p;
vtil = parin.vtil;
tolc = parin.tolc;
tolr = parin.tolr;
maxr = parin.maxr;

%% 
[a,b,V] = lanczosing(Schol,k+p,vtil);
s_old = ones(k,1);
T = diag(a) + diag(b(2:end-1),1) + diag(b(2:end-1),-1);
[Y,s_new] = eig(T);s_new = diag(s_new);
[s_new,ides] = sort(s_new,'descend');
Y = Y(:,ides);
nrest = 0;
sdiff = zeros(1,maxr); % record the "convergence history"
while norm(s_new(1:k) - s_old) > tolr && nrest < maxr
    sdiff(nrest + 1) = norm(s_new(1:k) - s_old);
    s_old = s_new(1:k);
    g = coeff(Y(:,1:k),k,p,tolc);
    vtil = V(:,1:k+p)*(Y(:,1:k)*g');
    [a,b,V] = lanczosing(Schol,k+p,vtil);
    T = diag(a) + diag(b(2:end-1),1) + diag(b(2:end-1),-1);
    [Y,s_new] = eig(T);s_new = diag(s_new);
    [s_new,ides] = sort(s_new,'descend');
    Y = Y(:,ides);
    nrest = nrest + 1;
end
sT = size(T,1);
ei = sqrt(s_new(1:k));
ei = 1./ei;
VQ1 = V(:,1:sT)*Y(:,1:k);
U = horzcat(VQ1,-Jmul(Schol'\(Schol\VQ1)));
if nrest == maxr
    warning('Restarting seems to fail')
    mess = 'Restarting seems to fail';
else
    mess = 'Restarting successfully';
end
sdiff = sdiff(1:nrest);
parout = struct;
parout.mess = mess;
parout.sdiff = sdiff;
parout.nrest = nrest;
end

