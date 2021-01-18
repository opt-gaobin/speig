function [S,D,Q,G] = williasondiag(A)
% This function computes a Williamson's diagonal form of spd A 
%  S^TAS = diag(D,D)
% where D = diag(d_1,...,d_n), d_j<=d_j+1 and S is symplectic
% The decomposition is based based on real Schur decomposition of A
% followed by rearrangment and multiplication with perfect shuffle permutation matrix.
% INPUT
% A: must be spd of even dimension
% OUTPUT
% S: the symplectic diagonalizing matrix 
% D: the symplectic eigenvalues of A
% Q: such that Q^TAQ = [0 D;-D 0]
% G: such that G^TAG = diagblk([0 dj,-dj 0])
% AUTHOR
% NTSON, UC Louvain, 28/01/2020
% =============REFERENCES ======================
% K.R. Parthasarathy, The symplectry group of Gaussian states in L^2(R^n),
% in A.N. Shiryaev et al. (eds.), Prokhorov and Contemporary Probability
% Theory, Springer Preceedings in Mathematics and Statistics 33, 
% DOI 10.1007/978-3-642-33549-5 21, © Springer-Verlag Berlin Heidelberg 2013
A = full(A);
if mod(size(A,1),2) ~= 0
    error('A must be even dimensional')
else
    [~,sh] = chol(A);
    if sh > 0
        error('A is not spd')
    else
        sqA = sqrtm(A);
        n = size(A,1)/2;
        JsqA = [sqA(n+1:2*n,1:n) sqA(n+1:2*n,1+n:2*n);...
            -sqA(1:n,1:n) -sqA(1:n,1+n:2*n)];
        B = sqA'*JsqA; B = 0.5*(B-B');
        [U,D]= schur(B);
        eshit = zeros(n,1);
        esign = zeros(1,n);
        for j = 1:n
            if D(2*j-1,2*j) < 0
                esign(j) = -1;
            end
            eshit(j) = abs(D(2*j-1,2*j));
        end
        [~,Ind] = sort(eshit,'ascend');
        esign = esign(Ind);
        G = zeros(2*n,2*n);
        % arranging block in increasing order
        for j = 1:n
            G(:,2*j-1:2*j) = U(:,2*Ind(j)-1:2*Ind(j));
        end
        % twisting inside blocks those have negative element above main diagonal 
        for j = 1:n
            if esign(j) < 0
                sh = G(:,2*j-1); G(:,2*j-1) = G(:,2*j); G(:,2*j)=sh;
            end
        end
        Q = zeros(2*n,2*n);
        for j = 1:n
            Q(:,j) = G(:,2*j-1);
            Q(:,n+j) = G(:,2*j);
        end
        D = sort(eshit,'ascend');
        sqD = 1./sqrt(D);
        S = (-JsqA*Q)*[zeros(n,n) diag(sqD);-diag(sqD) zeros(n,n)];
    end
end
end
