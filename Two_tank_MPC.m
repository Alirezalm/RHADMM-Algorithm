clc, clear, close all

%% time structure

t0 = 0;   tf = 500;        dt = 5;
time = t0: dt: tf - dt; Nt = numel(time);

%% System parameters
Aa = 0.06;
a1 = 6.7371e-4;
a3 = 4.0423e-4;
gama = 0.4;
g = 9.8;
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
KK = 0;
 % predection horizon

 alpha = 0.9;
    for    N1 = 5
        KK = KK + 1;
        for i = 1: N    % cost of MPC
            
            Q{i} = 20 * [1 0; 0 0.001];
            
            P{i} = Q{i};
            
            R{i} = 1 * eye(nu);
        end
        %% Local Constraints
        
        ubx = [1.36 1.30]';
        lbx = -ubx;
        Aineqx = [eye(nx);-eye(nx)];
        bineqx = [ubx; -lbx];
        ubu = 4;
        lbu = -ubu;
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
        
        b = 0.5 * ones(size(Atilda',1),1); % vector of Coupling Constraint
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
        
%         for i = 1:N
            x(:,1) = [1 0.26]';
            x(:,2) = [0.7 0.18]';
            x(:,3) = [0 0]';
            x(:,4) = [0.1 0.26]';
            x(:,5) = [0.32 0.26]';
            x(:,6) = [0.32 0.26]';
            %     x(:,i) =[unifrnd(lbx(1),ubx(1),1,1) unifrnd(lbx(2),ubx(2),1,1)]';
%         end
        %% Start Simulation
        
        
        ref = [2  0]';
        
        ref = 0.6 * (time <= 180) + 1.2 * (time <= 320 & time >= 180.1) + (0.4 *(time >= 320.1))  ;
        
        % ref = ref * 0;
        
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
            %     clc
            tic
            [xopt, fval, Data , theta  ] = RHCADMM(QQ,c,Aeq,beq,Aineq,bineq,M,b,A1,B1,alpha,N,m,x,options);
            GG(k) = toc;
            % GG(k)
            MM(k) = Data.Maxiter;
            % MM
            objective = 0;
            for i = 1:N
                objective = objective + 0.5 * theta(:,i)' * QQ{i} *  theta(:,i) ;    % Computing Optimal value
            end
            for i = 1: N
                
                u(:,i) = theta((N1+1)*nx + 1,i);
            end
            
            
            for i = 1: N
                x(:,i) = A{i} * x(:,i) + B{i} * u(:,i);
                
            end
            %      M{1}*theta(:,1)+M{2}*theta(:,2)+M{3}*theta(:,3)+M{4}*theta(:,4)+M{5}*theta(:,5)+M{6}*theta(:,6)
            RHO(k) = Data.finalrho;
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
            U8 (k) = U1 (k) + U2 (k) + U3 (k) + U4 (k);
            Objective_MPC(k) = objective;
            TIME(k) = toc;
            PP(k) = Data.Maxiter;
        end
        
        worstcase = max(GG);
        worstiter= max(MM);
        
        disp(['Horizon = ' num2str(N1) ' $\alpha$ = ' num2str(alpha) ' worst time = ' num2str(worstcase) ' worst iter = ' num2str(worstiter)])
end
        
        
N = 5:5:40;

alpha1 = 1000 * [0.511 0.520 0.533 0.583 0.700 0.876 1.037 1.183];

alpha2 = 1000 * [0.414 0.442 0.457 0.513 0.622 0.748 0.879 1.028];

alpha3 = 1000 * [0.335 0.382 0.421 0.502 0.596 0.726 0.821 1.080];

alpha4 = 1000 * [0.307 0.368 0.390 0.510 0.518 0.642 0.740 0.872];

alpha5 = 1000 * [0.286 0.328 0.353 0.394 0.464 0.566 0.659 0.853];


figure

plot(N, alpha1,N, alpha2,'-.*',N, alpha3, '->',N, alpha4,'-*',N, alpha5,'-o','linewidth',1.25)

legend('HADMM     \alpha = 0.5','RH-ADMM \alpha = 0.6','RH-ADMM \alpha = 0.7','RH-ADMM \alpha = 0.8','RH-ADMM \alpha = 0.9')

xlabel('Prediction Horizon (T)')

ylabel('Worst execution time (ms)')
%



figure
subplot(2,1,1)
plot(time,X1(1,:) ,'r --','linewidth',1.25)
hold on
plot(time,X2(1,:),'g :','linewidth',1.75)
hold on
plot(time,X3(1,:),'m -.','linewidth',1.25)
hold on
plot(time,X4(1,:),'c --','linewidth',1.25)
% hold on
% plot(time,X5(1,:),'m -.','linewidth',1.25)
% hold on
% plot(time,X6(1,:),'c --','linewidth',1.25)
hold on
plot(time,ref,'b --','linewidth',1.25)
% plot(time,0.68 * ones(Nt,1),'linewidth',1)
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
plot(time,X1(2,:) ,'r --','linewidth',1.25)
hold on
plot(time,X2(2,:),'g :','linewidth',1.75)
hold on
plot(time,X3(2,:),'m -.','linewidth',1.25)
hold on
plot(time,X4(2,:),'c --','linewidth',1.25)
hold on
% plot(time,X5(2,:),'m -.','linewidth',1.25)
% hold on
% plot(time,X6(2,:),'c --','linewidth',1.25)
hold on
plot(time,ref,'b--','linewidth',1.25)
% plot(time,0.65 * ones(Nt,1),'linewidth',1)
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
plot(time,U1 + 2,'r --','linewidth',1.25)
hold on
plot(time,U2+ 2,'g :','linewidth',1.75)
hold on
plot(time,U3+ 2,'m -.','linewidth',1.25)
hold on
plot(time,U4+ 2,'c --','linewidth',1.25)
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
plot(time,U8 + 4 * 2,'r --','linewidth',1.25)   % N = 4 ----> Nnmber of conected tank
hold on
plot(time,8.5*ones(1,Nt),'g .-','linewidth',1.25)
hold on
legend('$\sum_{i = 1}^{N} q^i$','total input flow')
xlabel('time [sec]')
ylabel('Coupling Constraint')

figure
plot(X1(1,:),X1(2,:),'b -.')
hold on
plot(X2(1,:),X2(2,:),'g --')
hold on
plot(X3(1,:),X3(2,:),'r -.')
hold on
plot(X4(1,:),X4(2,:),'c --')
hold on
plot(ref,ref,'m :','linewidth',1.5)
hold on
plot(X1(1,1),X1(2,1),'b *')
hold on
plot(X1(1,100),X1(2,100),'b o','markersize',11,'markerfacecolor','b')
hold on
plot(X2(1,1),X2(2,1),'g *')
hold on
plot(X2(1,100),X2(2,100),'g o','markersize',9,'markerfacecolor','g')
hold on
plot(X3(1,1),X3(2,1),'r *')
hold on
plot(X3(1,100),X3(2,100),'r o','markersize',7,'markerfacecolor','r')
hold on
plot(X4(1,1),X4(2,1),'c *')
hold on
plot(X4(1,100),X4(2,100),'c o','markersize',5,'markerfacecolor','c')
hold on
% plot(0.6,0.6,'m *')
% hold on
plot(0.4,0.4,'m o','markersize',3,'markerfacecolor','m')
% hold on
% plot(1.2,1.2,'m s','markersize',5,'markerfacecolor','m')
xlabel('x_1')
ylabel('x_2')
legend('tank 1', 'tank 2','tank 3', 'tank 4', 'Ref')
