% Testing codes for spd matrix with known symplectic eigenvalues
% References
% [1] P. Amodio, On the computation of few eigenvalues of positive definite
% Hamiltonian matrices, Future Generation Computer Systems 22(2006) 403-411
% [2] SAGS (2021) Symplectic eigenvalue problem via trace minimization and
% Riemannian optimization
% Author(s)
% NT Son, TNUS, 2020-12-17.
% NT Son, 2020-12-19. Used symplectic but not orthogonal matrix M
%% creating spd matrix
list_n = 100:100:500;
ln = length(list_n);
err = zeros(3,ln); % record the error
rest = zeros(3,ln); % record the residual;
eig_cell = cell(3,ln);
out_spopt = cell(1,ln);
out_splLanz = cell(1,ln);
k = 5;
% parameters for spopt
opts = struct;
opts.record = 0;
opts.mxitr  = 2000;
opts.xtol = 1e-11;
opts.ftol = 1e-11;
opts.gtol = 1e-9;
opts.maxtau = 1;
% parameters for rest_lanczos
par = struct;
par.p = k; % parameter for restarting
par.tolr = 1e-8; % tolerance for restarting
par.maxr = 2*k;  % the maximal number of restartings
par.tolc = 1e-7; % tolerance for comupting the coefficients
format long
disp('        eigs             symplLanczos          Riemannian')
for i = 1:ln
    n = list_n(i);
    % parameters for symplectic Gauss transformation
    Gauss.c = 1.2;
    Gauss.m = round(n/5);
    Gauss.d = -sqrt(Gauss.m);
    DD = diag([1:n 1:n]);
    rng default
    Q1 = randn(n,n); Q2 = randn(n,n);
    [U,~,~] = svd(Q1+sqrt(-1)*Q2); 
    M = [real(U) -imag(U);imag(U) real(U)];
    d1 = [ones(1,Gauss.m-2) Gauss.c Gauss.c ones(1,n-Gauss.m)];
    d2 = [ones(1,Gauss.m-2) 1/Gauss.c 1/Gauss.c ones(1,n-Gauss.m)];
    L2 = zeros(n,n); L2(Gauss.m,Gauss.m-1) = Gauss.d;
    L2(Gauss.m-1,Gauss.m) = Gauss.d;
    Lmcd = [diag(d1) L2;zeros(n,n) diag(d2)];%sympl. Gauss transformation type I
    M = M*Lmcd;
    M = 0.5*((M*DD)*M' + M*(DD*M'));
    M = 0.5*(M+M');
    %% eigs
    H = sparse([M(n+1:2*n,:);-M(1:n,:)]);
    [Veigs,eHs] = eigs(H,2*k,'SM');
    HVeigs = H*Veigs;
    rest(1,i) = norm(HVeigs - Veigs*eHs,'fro')/...
        norm(HVeigs,'fro');
    eHs = diag(eHs);
    eHs = eHs(2:2:2*k);
    eig_cell{1,i} = -sqrt(-1)*eHs;
    err(1,i)= norm(eig_cell{1,i}-(1:k)',1);
    %% restart Lanczos
    [eig_cell{2,i},U,out_splLanz{i}] = rest_lanczos(M,k,par);
    err(2,i)= norm(eig_cell{2,i}-(1:k)',1);
    [Mchol,chM] = chol(M,'lower');
    JMU = Jmul(Mchol'\(Mchol\U));
    rest(2,i) = norm(JMU-U*Jmul(diag([ones(k,1);1./...
        eig_cell{2,i}.^2])),'fro')/norm(JMU,'fro');
    %% optimizing
    X0 = zeros(2*n,2*k); X0(1:k,1:k) = eye(k);
    X0(n+1:n+k,k+1:end) = eye(k);
    [Xmin, out_spopt{i}]= spopt(X0, @eigvalcost, opts, M);
    MXmin = M*Xmin;
    XMXmin = Xmin'*MXmin;  XMXmin = 0.5*(XMXmin+XMXmin');
    [Sr,eig_cell{3,i}] = williasondiag(XMXmin);
    err(3,i) = norm(eig_cell{3,i}-(1:k)',1);
    Xmin = Xmin*Sr;MXmin = MXmin*Sr;
    JXmin = [Xmin(n+1:2*n,:);-Xmin(1:n,:)];
    rest(3,i) = norm(MXmin - JXmin*[zeros(k,k) -diag(eig_cell{3,i});...
        diag(eig_cell{3,i}) zeros(k,k)],'fro')/...
    norm(MXmin,'fro');
disp(['matrix size ',num2str(2*n)])
disp(horzcat(real(eig_cell{1,i}), eig_cell{2,i}, eig_cell{3,i}))
end
%
figure(1)
semilogy(list_n,err(1,:),'b-->',...
    list_n,err(2,:),'r-.o',...
    list_n,err(3,:),'k:s',...
    'LineWidth',1.5)
legend('eigs','symplLanczos','Riemannian','Location','east')
xlabel('problem size n')
ylabel('errors')
set(gcf,'color','w')



figure(2)
semilogy(list_n,rest(1,:),'b-->',...
    list_n,rest(2,:),'r-.o',...
    list_n,rest(3,:),'k:s',...
    'LineWidth',1.5)
legend('eigs','symplLanczos','Riemannian','Location','southeast')
xlabel('problem size n')
ylabel('residuals')
set(gcf,'color','w')
%}
format short
%% ============= functions ==========================
function [F,G] = eigvalcost(X,M)
MX = M*X;
F = X(:)'*MX(:);
G = 2*MX;
end
