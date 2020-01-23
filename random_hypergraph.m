function [H,A,B] = random_hypergraph(M,N,R,overlap)
% RANDOM_HYPERGRAPH  Build a random geometric graph.
%   [H,A,B] = RANDOM_HYPERGRAPH(M,N,R,overlap) generates a network with M
%   fusion centres, N nodes and maximum distance to connect to a FC equal
%   to R. overlap is true if the FC are nodes, false if they are separate
%   agents. The routine outputs a cell array with elements representing the
%   hyper-edges, and matrices A and B that appear in the constraints.


%% generate network
% generate nodes-FCs incidence matrix
I = zeros(N,M);

flag = true;
while(flag)
    
    % nodes' positions
    N_pos = rand(N,1);
    
    % position of FCs
    if overlap
        % positions of nodes chosen to be FCs
        F_pos = N_pos(randsample(N,M));
    else
        % random FCs positions
        F_pos = rand(M,1);
    end
    
    % compute distances between nodes and FCs
    for n = 1:N
        
        % compute distance from each FC
        I(n,:) = abs(N_pos(n) - F_pos) <= R;
        
    end
    
    % if each node is not connected to at least a FC, generate new N_pos
    flag = ~all(sum(I,2));
    
end

%% A and B matrices
% find indices of node-FC connections
[N_idx,F_idx] = find(I);

% number of constraints
T = length(N_idx);

A = zeros(T,N);
B = zeros(T,M);

for t = 1:T
    
    A(t,N_idx(t)) = -1;
    B(t,F_idx(t)) = 1;
    
end

%% hyper-edge cell array
H = cell(M,1);

for m = 1:M
    
    H{m} = find(I(:,m));
    
end


end

