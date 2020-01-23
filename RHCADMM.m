function [xopt,fval, Data , theta  ] = RHCADMM(Q,c,Aeq,beq,Aineq,bineq,M,b,A,B,alpha,N,m,state,options)

nx = size(Q{1} , 1);             % dimension of primal decision variable
n = size(M{1} , 1);              % dimension of the dual variable lambda(consensus problem)\\



T = size(A,1)/n;                % number of constraints (block-wise)

%% Simulation Parameters


MaxIter = options.MaxIter;

tol = options.tolerance;

eps = options.hypergraph;

if strcmpi(options.Display,'OFF')
    
    Quiet = 0;
    
elseif strcmpi(options.Display,'ON')
     Quiet = 1;
end







%% Fix Penalty parameter selection

rho = 1;   % without adaptive penalty




%% Offline computational parts
E = B' * B;
%% OFF-line Part of Computations
inv_E = E^-1;

H = B * inv_E * B';

z = 1*ones(T * n , 1);
% x = zeros(N * n , 1);
y = 1*ones(n * m , 1);
Ed = diag(B' * B);    % Computing e_{j}

Dd = diag(A' * A);    % Computing d_{j}

e = zeros(1,m);
d = zeros(1,N);

for j = 1: m
    e(j) = Ed((j-1) * n + 1);
end

for i = 1: N
    d(i) = Dd((i-1) * n + 1);
end
    if Quiet
%         clc
%     disp(['iter = ' num2str(k) '  Pri_Res = ' num2str(r) '  Dual_Res = ' num2str(s) '  Objective = ' num2str(f)])
   fprintf('%3s\t\t%4.6s\t\t%4.6s\t\t%4.6s\n','Iter','P_Res','D_RES',' Obj');
    end

%% ON-line Part of Computations (in Parallel)

for k = 1: MaxIter
    
      Data.steplength(k) = rho;
    
   [x , z , theta, fval] = xupdate(Q,c,Aeq,beq,Aineq,bineq,M,b,rho,alpha,z,A,B,eps,n,N,d,m,H,inv_E,state );  % xupdate function computes the x-updates and z-updates parts
    yold = y;
    y = ((1/rho)*max(E^-1*B' * z,0));                               % Computing y
    s = norm(rho * A' * B * (y - yold));
    r = norm (A * x(:) + B * y);       
    if alpha ~= 0.5 % Computing Resedual of the equality constraint
    rho = linesearch(rho,x,yold,z,A,B,inv_E);
    end
    f = 0;
    for i = 1:N
        f = f + 0.5 * theta(:,i)' * Q{i} *  theta(:,i) + c{i} * theta(:,i);    % Computing Optimal value
    end
    %% Saving Results
    
    Data.PrimalResedual(k) = r;
    Data.DualResedual(k) = s;
    Data.Maxiter = k;
    Data.Cost(k) = sum(fval);
   
    
        %% Display the results
    
    if Quiet
%         clc
%     disp(['iter = ' num2str(k) '  Pri_Res = ' num2str(r) '  Dual_Res = ' num2str(s) '  Objective = ' num2str(f)])
    fprintf('%3d\t\t%4.6f\t\t%3.3f\t\t%3.3f\n', k, ...
            r,s, ...
            fval);
    end
    %% Stopping ceriterion 
    if r <= tol && s <= tol
        Data.Maxiter = k;
        xopt = x(:,1);
        fval = f;
        Data.finalrho = rho;
%         disp('The optimization problem has been solved') 
        break
    elseif k == MaxIter
        Data.Maxiter = k;
        xopt = x(1:n);
        fval = f;
        disp('The algorithm is terminated due to the number of iterations') 
        break
        
    end
end

end

