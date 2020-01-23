function  [x , z, theta,fval] = xupdate(Q,c,Aeq,beq,Aineq,bineq,M,b,rho,alpha,z,A,B,eps,n,N,d,m,H,inv_E,state)
opt = optimoptions('quadprog', 'Display', 'off', 'algorithm','interior-point-convex');
n1 = zeros(1,m);
for i = 1:m
    n1(i) = size(eps{i},2);
end
beta = ((z - max(2 * H * z , 0))' * A)';
%% Inner QP solve
nx = size(Q{1},1);

% for i = 1: N
%     theta(:,i) = quadprog(0.5*((Q{i} + (1/(rho*d(i)))*M{i}'*M{i})+(Q{i} + (1/(rho*d(i)))*M{i}'*M{i})'), -(1/(rho*d(i)))*...
%         ((b/N)'*M{i} - beta(n*(i - 1) + 1:i*n)'*M{i} - c{i}),[],[],Aeq{i},beq{i} * state(:,i),[],[],[],opt);
%
% % end
% for i = 1: N
%     HH{i} = Q{i} + (1/rho) * M{i}' * (d(i) * eye(n))^-1 * M{i};
%     HH{i} = 0.5 * ( HH{i} +  HH{i}');
%     grad{i} = (-c{i} + (1/rho) * ( (b/N)' * (d(i) * eye(n))^-1  - beta(n*(i - 1) + 1:i*n)' * (d(i) * eye(n))^-1 )* M{i}   ) ;
% %     [theta(:,i),fval(i)] = quadprog(HH{i},-grad{i},Aineq{i},bineq{i},Aeq{i}, beq{i} * state(:,i),[],[],[],opt);
% %  theta(:,i) = thetaa;
% end
% tic
% parfor i = 1: N
% %     HH{i} = Q{i} + (1/rho) * M{i}' * (d(i) * eye(n))^-1 * M{i};
% %     HH{i} = 0.5 * ( HH{i} +  HH{i}');
% %     grad{i} = (-c{i} + (1/rho) * ( (b/N)' * (d(i) * eye(n))^-1  - beta(n*(i - 1) + 1:i*n)' * (d(i) * eye(n))^-1 )* M{i}   ) ;
%     [theta(:,i),fval(i)] = quadprog(HH{i},-grad{i},Aineq{i},bineq{i},Aeq{i}, beq{i} * state(:,i),[],[],[],opt);
% %  theta(:,i) = thetaa;
% end
% toc
% tic
for i = 1: N
    HH{i} = Q{i} + (1/rho) * M{i}' * (d(i) * eye(n))^-1 * M{i} + 0.1 * eye(size(Q{i},1));
    HH{i} = 0.5 * ( HH{i} +  HH{i}');
    grad{i} = (-c{i} + (1/rho) * ( (b/N)' * (d(i) * eye(n))^-1  - beta(n*(i - 1) + 1:i*n)' * (d(i) * eye(n))^-1 )* M{i}   ) ;
    [thetaa,fval] = quadprog(HH{i},-grad{i},Aineq{i},bineq{i},Aeq{i}, beq{i} * state(:,i),[],[],[],opt);
 theta(:,i) = thetaa;
end
% toc
for i = 1:N
    x(:,i) = (1/rho) * (d(i) * eye(n))^-1 * (M{i} * theta(:,i) - (b/N) + beta(n*(i - 1) + 1:i*n));
end


xx = x(:);

z = z - max(2 * alpha * B * inv_E * B' * z , 0) - 2 * alpha * rho * A * xx;


z = z(:);
end

