function  rho = linesearch(rho,x,y,z,A,B,inv_E)

tau = 3;
miu = 10;
yold = y;
y = (1/rho) * max(inv_E*B' * z,0);                               % Computing y



r =  A * x(:) + B * y;

S = rho * A' * B * (y - yold);

% if norm(S) ~= 0
    if norm(r) > miu * norm(S)
        rho = rho*(1 + tau);
        
    elseif norm(S)> miu * norm(r)
        
        rho = rho/(1 + tau);
        
    end
% end


end



