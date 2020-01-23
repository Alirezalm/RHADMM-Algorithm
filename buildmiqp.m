function [H,G,Ibar,Qbar,DD,dd] = buildmiqp(N,Q,P,R,A,B,Aineqx,bineqx,Ainequ,binequ)
nx = size(A,1); nw = size(B,2);  nix = size(Aineqx,1); niu = size(Ainequ,1); 
for i = 1 : N
    Qbar1(nx*i-nx+1:nx*i ,nx*i-nx+1:nx*i ) = Q ;
end

for i = 1 : N
    Qbar1(nx*i-nx+1:nx*i ,nx*i-nx+1:nx*i ) = Q ;
end

Qbar = blkdiag(Qbar1,P);


for i = 1 : N
    Rbar(nw*i-nw+1:nw*i ,nw*i-nw+1:nw*i ) = R ;
end

H = 2 * blkdiag(Qbar,Rbar);


for i = 1 : N
    Abar(nx*i-nx+1:nx*i ,nx*i-nx+1:nx*i ) = A ;
end
Abar = [zeros(nx,N*nx);Abar];
Abar = [Abar zeros((N+1)*nx,nx)];

for i = 1 : N
    Bbar(nx*i-nx+1:nx*i ,nw*i-nw+1:nw*i ) = B ;
end
Bbar = [zeros(nx,N*nw); Bbar];

Ibar = eye(nx);
Ibar = [Ibar;zeros((N)* nx , nx)];



G = [eye(size(Abar,1))-Abar -Bbar];


for i = 1 : N + 1
    DDx(nix*i-nix+1:nix*i ,nx*i-nx+1:nx*i ) = Aineqx ;
end

for i = 1 : N + 1
    ddx(nix*i-nix+1:nix*i ,: ) = bineqx ;
end

for i = 1 : N
    DDu(niu*i-niu+1:niu*i ,nw*i-nw+1:nw*i ) = Ainequ ;
end

for i = 1 : N
    ddu(niu*i-niu+1:niu*i ,: ) = binequ ;
end

DD = blkdiag(DDx,DDu);
dd = [ddx;ddu];

end

