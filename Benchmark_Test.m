%% Benchmark Test for the distributed QP problem of the following form

clc, clear, close all

%%  Structure of the hypergraph

N = 6;                        % Number of Nodes

m = 5;                        % Number of Fusion Centers

% Case I: N = 6, m = 5
eps{1} = [1 4];

eps{2} = [1 3];

eps{3} = [2 3];

eps{4} = [2 5];

eps{5} = [5 6];


%% Generate problem suitable Positive Definite Matrices Q_{i}s and vectors c_{i}s and constraint matrices A and B

n = 500;                        % Dimension of x_{i}
mm = 5;
for i = 1 : N
    
    H = rand(n);
    
    H = 0.5*(H + H');
    
    [V, ~] = eig(H);
    
    Q{i}= V*diag(1+rand(n,1))*V';
    c{i} = 1*randn(n,1);
    M{i} = eye(mm,n);
 
end
b = 1* ones(mm,1);
% Case one
A = [-eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) -eye(mm) zeros(mm) zeros(mm);
    -eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) -eye(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) -eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) -eye(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) -eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) zeros(mm) -eye(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) zeros(mm) -eye(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm) -eye(mm)];
B = [ eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
      eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
      zeros(mm) eye(mm) zeros(mm) zeros(mm) zeros(mm);
      zeros(mm) eye(mm) zeros(mm) zeros(mm) zeros(mm);
      zeros(mm) zeros(mm) eye(mm) zeros(mm) zeros(mm);
      zeros(mm) zeros(mm) eye(mm) zeros(mm) zeros(mm);
      zeros(mm) zeros(mm) zeros(mm) eye(mm) zeros(mm);
      zeros(mm) zeros(mm) zeros(mm) eye(mm) zeros(mm);
      zeros(mm) zeros(mm) zeros(mm) zeros(mm) eye(mm);
      zeros(mm) zeros(mm) zeros(mm) zeros(mm) eye(mm)];


maxiter = 30;

tol = 1e-3;

options.MaxIter = maxiter;

options.tolerance = tol;

options.hypergraph = eps;

options.Display = 'on';

%% The Relaxed Hybrid Consensus ADMM (RHCADMM) solver solves the problem  

load ali Q c
 alpha = 0.9;
[x,~,~,~,lambda] = quadprog(blkdiag(Q{1},Q{2},Q{3},Q{4},Q{5},Q{6}),[c{1}',c{2}',c{3}',c{4}',c{5}',c{6}'],[M{1},M{2},M{3},M{4},M{5},M{6}],b,[],[],0.1 * ones(N*n,1));


tic
[xopt ,Data,theta] = RHCADMM(Q,c,M,b,A,B,alpha,N,m,options);
toc
disp(['alpha = ' num2str(alpha) '  maxiter = ' num2str(Data.Maxiter)])





subplot(3,1,1)
plot(1:Data.Maxiter, Data.PrimalResedual)
hold on
plot(1:Data.Maxiter, Data.DualResedual)
hold on
ylabel('Resedual')
xlabel('iterations')
legend('Primal Resedual', 'Dual Resedual')

subplot(3,1,2)
plot(1:Data.Maxiter, Data.Cost)
hold on

ylabel('Optiaml Value')
xlabel('iterations')

subplot(3,1,3)
plot(1:Data.Maxiter,Data.steplength)
hold on

ylabel('\rho')
xlabel('iterations')
% end
% theta
% % fval1
% % fval