function JA = Jmul(A)
% This function performs multiplication of J and A where J is the Poisson
% matrix
% Author(s)
% NT Son, UCLouvain, 2020-02-12
rA = size(A,1)/2;
JA = [A(rA+1:2*rA,:);-A(1:rA,:)];
end

