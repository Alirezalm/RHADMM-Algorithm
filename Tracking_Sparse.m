clc, clear, close all

%% time structure

t0 = 0;         tf = 5;        dt = 0.1;
time = t0: dt: tf - dt; Nt = numel(time);


%% System Structure
N = 6;
m = 5;
% A{1} = [0 1;-8 -3];   B{1} = [0 1]';
% A{2} = [0 1;-1 -2];   B{2} = [0 1]';
% A{3} = [0 1;-1 -1];   B{3} = [0 1]';
% A{4} = [0 1;-1 -5];   B{4} = [0 1]';
% A{5} = [0 1;-1 -7];   B{5} = [0 1]';
% A{6} = [0 1;-1 -8];   B{6} = [0 1]';
for i = 1:N
%     A{i} = [0 1;-1 -2];   B{i} = [0 1]';
%     
%     [A{i}, B{i}] = c2d(A{i},B{i},dt);
    A{i} = [0.875 0.1250; 0.1250 0.8047];   B{i} = [0.3 3]';
end


nx = size(B{1},1);

nu = size(B{1},2);
%% MPC Parameters
N1 = 5;

for i = 1: N
    
    Q{i} = 20 * rand *[1 0; 0 1];
    
    P{i} = Q{i};
    
    R{i} = 5 * rand * eye(nu);
end
%% Local Constraints

ubx = [1 0.64]';
lbx = - ubx;
Aineqx = [eye(nx);-eye(nx)];
bineqx = [ubx; -lbx];
ubu = 0.3;
lbu = - ubu;
Ainequ = [eye(nu);-eye(nu)];
binequ = [ubu; -lbu];

%% DMPC Matrices
for i = 1: N
    
    [H,G,Ibar,Qbar,DD,dd] = buildmiqp(N1,Q{i},P{i},R{i},A{i},B{i},Aineqx,bineqx,Ainequ,binequ);
    
    QQ{i} = H;
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

b = 0.8* ones(size(Atilda',1),1);
mm = size(Atilda',1);
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
%% Initial Conditions

for i = 1:N
    x(:,i) =[unifrnd(lbx(1),ubx(1),1,1) unifrnd(lbx(2),ubx(2),1,1)]';
end
%% Start Simulation


ref = [2  0]';

ref = 0.8* (time <= 3) + 0.3 * (time <= 6 & time >= 3.1)+(-0.8*(time >= 6.1))  ;

ref = ref * 0;

% h = waitbar(0,'Simulation Progress');

k = 0;
for t = t0 : dt: tf - dt
    
    k = k + 1;
    r1 = [];
    for i = 1: N1 + 1
        
        r = [ref(k)*1 0]';
        
        r1 = [r;r1];
    end
    
    X1(:,k) = x(:,1);
    X2(:,k) = x(:,2);
    X3(:,k) = x(:,3);
    X4(:,k) = x(:,4);
    X5(:,k) = x(:,5);
    X6(:,k) = x(:,6);
    
    for i = 1: N
        f = [-2 * r1' * QQbar{i} zeros(1, N1 * nu)];
        c{i} = f;
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
plot(time,X1(1,:) + 1,'--','linewidth',1)
hold on
plot(time,X2(1,:)+ 1,'-*','linewidth',1)
hold on
plot(time,X3(1,:)+ 1,'-.','linewidth',1)
hold on
plot(time,X4(1,:)+ 1,'-o','linewidth',1)
hold on
plot(time,ones(Nt,1),'linewidth',1)
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
plot(time,X1(2,:) + 0.64,'--','linewidth',1)
hold on
plot(time,X2(2,:) + 0.64,'-*','linewidth',1)
hold on
plot(time,X3(2,:)+ 0.64,'-.','linewidth',1)
hold on
plot(time,X4(2,:)+ 0.64,'-o','linewidth',1)
hold on
plot(time,0.64 * ones(Nt,1),'linewidth',1)
xlabel('time [sec]')
ylabel('h_2 [m]')
legend('tank 1', 'tank 2','tank 3', 'tank 4', 'Ref')
% plot(time,X5(2,:)+ 0.64)
% hold on
% plot(time,X6(2,:)+ 0.64)
% hold on
% plot(time,ones(Nt,1)*lbx(1),time,ones(Nt,1)*ubx(1))
% plot(time,ref)
figure
subplot(2,1,1)
plot(time,U1 + 0.3,'--','linewidth',1)
hold on
plot(time,U2+ 0.3,'-*','linewidth',1)
hold on
plot(time,U3+ 0.3,'-.','linewidth',1)
hold on
plot(time,U4+ 0.3,'-o','linewidth',1)
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
plot(time,U7 + N * 0.3,'--','linewidth',1)
hold on
plot(time,2*ones(1,Nt),'--','linewidth',1)
hold on
legend('$\sum_{i = 1}^{N} q^i$','total input flow')
xlabel('time [sec]')
ylabel('Coupling Constraint')
