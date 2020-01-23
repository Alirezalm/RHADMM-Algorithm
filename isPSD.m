function tf = isPSD(Q)

if issymmetric(Q)
    
    lambda = eig(Q);
    
    if lambda >= 0
        
        tf = true;
    else
        tf = false;
    end
else
    
    
%     warning('The Hessian matrix is not symmetric. Resetting Q = 0.5 * (Q + Q)')
    
    Q = 0.5 * (Q + Q');
    lambda = eig(Q);
    if lambda >= 0
        
        tf = true;
    else
        tf = false;
    end
end

