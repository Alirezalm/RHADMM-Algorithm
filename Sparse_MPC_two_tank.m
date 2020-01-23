clc, clear, close all

%% time structure

t0 = 0;         tf = 300;        dt = 5;
time = t0: dt: tf - dt; Nt = numel(time);

%% System parameters
Aa = 0.06;
a1 = 6.7371e-4;
a3 = 4.0423e-4;
gama = 0.4;
g = 9.81;
h1_0 = 0.68;
h3_0 = 0.65;
q_0 = 2;

tao_1 = (Aa / a1) * (sqrt(2*h1_0/g));
tao_3 = (Aa / a3) * (sqrt(2*h3_0/g));

%% System Structure
N = 6;
m = 5;

for i = 1:N
    A{i} = [-1/tao_1 1/tao_3;0 -1/tao_3];   B{i} = [gama/Aa (1-gama)/Aa]';
    
    [A{i}, B{i}] = c2d(A{i},B{i},dt);
end

nx = size(B{1},1);   % number of state

nu = size(B{1},2);  %  number of input
%% MPC Parameters
N1 = 5;  % predection horizon

for i = 1: N    % cost of MPC
    
    Q{i} = 20 * [1 0; 0 1];
    
    P{i} = Q{i};
    
    R{i} =  0.5 * eye(nu);
end
%% Local Constraints

ubx = [1.36 1.30]';
lbx = [-1.36 -1.30]';
Aineqx = [eye(nx);-eye(nx)];
bineqx = [ubx; -lbx];
ubu = 4;
lbu = -4;
Ainequ = [eye(nu);-eye(nu)];
binequ = [ubu; -lbu];

%% DMPC Matrices
for i = 1: N
    [H,G,Ibar,Qbar,DD,dd] = buildmiqp(N1,Q{i},P{i},R{i},A{i},B{i},Aineqx,bineqx,Ainequ,binequ);
    % PROBlEM of QP 
    QQ{i} = H;     % HESIAN MATRIX
    Aeq{i} = G;
    beq{i} = Ibar;
    QQbar{i} = Qbar;   
    Aineq{i} = DD;
    bineq{i} = dd;
end

%% Coupling Constraint

Atilda = zeros((N1+1)*nx,N1* nu);
Btilda = eye(N1*nu);
for i = 1: N
    M{i} =  [Atilda' Btilda];
end

b = 12 * ones(size(Atilda',1),1); % vector of Coupling Constraint
mm = size(Atilda',1);    %size x=n
% opt = optimoptions('quadprog','Display','off');
%% Graph Structure and Fusion Centers

A1 = [-eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) -eye(mm) zeros(mm) zeros(mm);
    -eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) -eye(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) -eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) -eye(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) -eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) zeros(mm) -eye(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) zeros(mm) -eye(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm) -eye(mm)];
B1 = [ eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    eye(mm) zeros(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) eye(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) eye(mm) zeros(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) eye(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) eye(mm) zeros(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) eye(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) eye(mm) zeros(mm);
    zeros(mm) zeros(mm) zeros(mm) zeros(mm) eye(mm);
    zeros(mm) zeros(mm) zeros(mm) zeros(mm) eye(mm)];
alpha = 0.9;

eps{1} = [1 4];

eps{2} = [1 3];

eps{3} = [2 3];

eps{4} = [2 5];

eps{5} = [5 6];

%% RHCADMM Solver options
maxiter = 1000;

tol = 1e-2;

options.MaxIter = maxiter;

options.tolerance = tol;

options.hypergraph = eps;

options.Display = 'off';
%% Initial Conditions of ststes system

for i = 1:N
    x(:,i) =[unifrnd(lbx(1),ubx(1),1,1) unifrnd(lbx(2),ubx(2),1,1)]';
%  x(:,i) = [0.32 1.26]';
end
%% Start Simulation

ref = [2  0]';

ref = 0.8* (time <= 3) + 0.3 * (time <= 6 & time >= 3.1)+(-0.8*(time >= 6.1))  ;

ref = ref * 0;

% h = waitbar(0,'Simulation Progress');

%%% simulation of MPC
k = 0;
for t = t0 : dt: tf - dt
    
    k = k + 1;
    r1 = [];
    for i = 1: N1 + 1
        
        r = [ref(k)*1 0]';
        
        r1 = [r;r1];
    end
    
    % STOREGE of ststes
    X1(:,k) = x(:,1);
    X2(:,k) = x(:,2);
    X3(:,k) = x(:,3);
    X4(:,k) = x(:,4);
    X5(:,k) = x(:,5);
    X6(:,k) = x(:,6);
    
    for i = 1: N
        f = [-2 * r1' * QQbar{i} zeros(1, N1 * nu)];
        c{i} = f;  % c of objective function of dual probem
    end
    
%     tic
    [xopt, fval, Data , theta  ] = RHCADMM(QQ,c,Aeq,beq,Aineq,bineq,M,b,A1,B1,alpha,N,m,x,options);
% toc
         objective = 0;
    for i = 1:N
        objective = objective + 0.5 * theta(:,i)' * QQ{i} *  theta(:,i) + c{i} * theta(:,i);    % Computing Optimal value
    end
    for i = 1: N

        u(:,i) = theta((N1+1)*nx + 1,i);
    end

    
    for i = 1: N
        x(:,i) = A{i} * x(:,i) + B{i} * u(:,i);
        
    end
%      M{1}*theta(:,1)+M{2}*theta(:,2)+M{3}*theta(:,3)+M{4}*theta(:,4)+M{5}*theta(:,5)+M{6}*theta(:,6)

        clc
        disp(['simulation progress = ' num2str((t+ dt) * 100 / tf) '%'])
%     waitbar(t / tf)
    %% Saving Data
    
    U1 (k)= u(:,1);
    U2 (k)= u(:,2);
    U3 (k)= u(:,3);
    U4 (k)= u(:,4);
    U5 (k)= u(:,5);
    U6 (k)= u(:,6);
    U7(k) = sum(u);
    Objective_MPC(k) = objective;

end


% close(h)
figure
subplot(2,1,1)
plot(time,X1(1,:) + 0.68,'g --','linewidth',1)
hold on
plot(time,X2(1,:)+ 0.68,'r ','linewidth',1)
hold on
plot(time,X3(1,:)+ 0.68,'k -.','linewidth',1)
hold on
plot(time,X4(1,:)+ 0.68,'c','linewidth',1)
hold on
plot(time,0.68 * ones(Nt,1),'b--','linewidth',0.5)
xlabel('time [sec]')
ylabel('h_1 [m]')
% plot(time,X5(1,:)+ 1)
% hold on
% plot(time,X6(1,:)+ 1)
% hold on
% plot(time,ones(Nt,1)*0)
% hold on
% plot(time,ref)
% plot(time,ref)
legend('tank 1', 'tank 2','tank 3', 'tank 4', 'Ref')
subplot(2,1,2)
plot(time,X1(2,:) + 0.65,'g --','linewidth',1)
hold on
plot(time,X2(2,:) + 0.65,'r','linewidth',1)
hold on
plot(time,X3(2,:)+ 0.65,'k -.','linewidth',1)
hold on
plot(time,X4(2,:)+ 0.65,'c','linewidth',1)
hold on
plot(time,0.65 * ones(Nt,1),'linewidth',1)
xlabel('time [sec]')
ylabel('h_3 [m]')
legend('tank 1', 'tank 2','tank 3', 'tank 4', 'Ref')
% plot(time,X5(2,:)+ 0.64)
% hold on
% plot(time,X6(2,:)+ 0.64)
% hold on
% plot(time,ones(Nt,1)*lbx(1),time,ones(Nt,1)*ubx(1))
% plot(time,ref)


figure
subplot(2,1,1)
plot(time,U1 + 2,'--','linewidth',1)
hold on
plot(time,U2+ 2,'-*','linewidth',1)
hold on
plot(time,U3+ 2,'-.','linewidth',1)
hold on
plot(time,U4+ 2,'-o','linewidth',1)
hold on

% stairs(time,U5+ 0.3)
% hold on
% stairs(time,U6+ 0.3)
% hold on
% plot(time,ones(Nt,1)*0.6,'-','linewidth',1)
legend('tank 1', 'tank 2','tank 3', 'tank 4')
xlabel('time [sec]')
ylabel('q(t)')
subplot(2,1,2)
plot(time,U7 + N * 2,'--','linewidth',1)
hold on
plot(time,20 * ones(1,Nt),'--','linewidth',1)
hold on
legend('$\sum_{i = 1}^{N} q^i$','total input flow')
xlabel('time [sec]')
ylabel('Coupling Constraint')
