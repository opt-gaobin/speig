function v = Sorthog(v,Schol,V,W,a)
%global Sv
% This function S-orthogonalize vector v w.r.t the columns of V and W
% bla bla
% INPUT
% V, W: two given parts of a symplectic matrix U = [V W]
% v: any vector of the same size as the number of rows of V
% a: [alp_1,...,alp_k]
% Reference
% P. Amodio, On the computation of few eigenvalues of positive definite
% Hamiltonian matrices, Future Generation Computer Systems 22(2006) 403-411
% Author(s)
% NT Son, UCLouvain, 2020-03-11
% NT Son, 2020-03-28, compute the smallest
% NT Son, 2020-12-17: replace S with its Cholesky factor
for j = 1:length(a)
    %r = (W(:,j)'*(S*v))/a(j);
    r = (W(:,j)'*(Schol'\(Schol\v)))/a(j);
    v = v - r*W(:,j);
    %r = V(:,j)'*(S*v);
    r = V(:,j)'*(Schol'\(Schol\v)); 
    v = v - r*V(:,j);
end
end

