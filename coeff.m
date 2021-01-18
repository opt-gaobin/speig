function g = coeff(Y,k,p,tolcoeff)
% This function computes the coefficient needed in restarted Lanczos method
% Reference
% P. Amodio, On the computation of few eigenvalues of positive definite
% Hamiltonian matrices, Future Generation Computer Systems 22(2006) 403-411
% Author(s)
% NT Son, UCLouvain, 2020-02-12
j = 1;
while norm(Y(p+j+1:p+k)) < tolcoeff && j <= k-1
    j = j + 1;
end
if j == k-1
    warning('coeff.m seems to fail')
end
Am = [ones(1,k-j+1);Y(p+j+1:p+k,j:k)];
g = [ones(1,j-1) [k+1-j zeros(1,k-j)]/Am'];
end

