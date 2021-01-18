function [a,b,V] = lanczosing(Schol,p,Vtil,ak,bk)
% This is the Lanczos function producing matrix S^{-1}-orthogoanl V and the 
% element of tridiagonal elements of T such that
% -H^2V = VT + beta_k+1v_k+1e_k^T
% where H = JS^{-1}, T  = diag(a) + diag(b,1) + diag(b,-1);
% INPUT
%   Schol: the lower Cholesky factor of spd matrix S such that H = JS
%   p: the number of Lanczos steps performed 
%   Vtil: the initial guess, which can be a vector if the process starts at
%   nothing and a matrix whose last column is the initial vector and the
%   left is the already part of matrix V computed. 
%   ak,bk: the entries of matrix T which are already computed
%   length(ak)+1 = length(bk)+ 1 = size(Vtil,2);
% OUTPUT
%   a, b: the diagonal and subdiagonal of matrix T
%   V: the matrix.
% Reference
% P. Amodio, On the computation of few eigenvalues of positive definite
% Hamiltonian matrices, Future Generation Computer Systems 22(2006) 403-411
% Author(s)
% NT Son, UCLouvain, 2020-03-11
% NT Son, 2020-3-13: modified from lanczos_basic to be more adaptive with input
% NT Son, 2020-03-28: compute the smallest
% NT Son, 2020-12-17: replace S with its Cholesky factor
%% 
% note H = JS
if nargin < 3
    error('function needs at least 3 inputs')
end
if nargin < 4; ak = []; bk = []; end
if size(Vtil,2) ~= length(ak)+1 || size(Vtil,2) ~= length(bk)+1
    error('check the compatibility of Vtil,ak,bk')
end
k = length(ak);
vtil = Vtil(:,end);
b = [bk zeros(1,p+1)];
a = [ak zeros(1,p)];
V = [Vtil(:,1:end-1) zeros(size(Schol,1),p+1)];
%W = [-Jmul(S*Vtil(:,1:end-1)) zeros(size(S,1),p)];
%b(k+1) = sqrt(vtil'*S*vtil);
W = [-Jmul(Schol'\(Schol\Vtil(:,1:end-1))) zeros(size(Schol,1),p)];
b(k+1) = sqrt(vtil'*(Schol'\(Schol\vtil)));
V(:,k+1) = vtil/b(k+1);
for j = 1:p
    %W(:,k+j) = -Jmul(S*V(:,k+j));
    %Sw = S*W(:,k+j);
    W(:,k+j) = -Jmul(Schol'\(Schol\V(:,k+j)));
    Sw = Schol'\(Schol\W(:,k+j));
    a(k+j) = W(:,k+j)'*Sw;
    if j == 1 && k == 0
            vtil = Jmul(Sw) - a(k+j)*V(:,k+j);
    else
        vtil = Jmul(Sw) - a(k+j)*V(:,k+j)- b(k+j)*V(:,k+j-1);
        vtil = Sorthog(vtil,Schol,V(:,1:k+j),W(:,1:k+j),a(1:k+j));
    end
    %b(k+j+1) = sqrt(vtil'*S*vtil);
    b(k+j+1) = sqrt(vtil'*(Schol'\(Schol\vtil)));
    V(:,k+j+1) = vtil/b(k+j+1);
end
end

