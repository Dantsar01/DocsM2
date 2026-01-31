function A = LapNeu2D(Nx)

I = speye(Nx,Nx);
u = ones(Nx-1,1);
H = 4*I - diag(u,1) - diag(u,-1);
H(1,1) = 1;
H(Nx,Nx) = 1;
H = sparse(H);

I2 = I;
I2(1,1) = 0;
I2(Nx,Nx) = 0;

A = sparse(Nx^2,Nx^2);
v= 1:Nx;
for i = 2 : Nx-1
    A((i-1)*Nx+v,(i-1)*Nx+v) = H;
    A((i-1)*Nx+v,(i)*Nx+v) = -I2;
    A((i-1)*Nx+v,(i-2)*Nx+v) = -I2;
end


% Condition de Neumann au bord
A(v,v) = I;
A(v,Nx+v) = -I;

A(Nx^2 - Nx  +v,Nx^2 - Nx +v) = I;
A(Nx^2 - Nx  +v,Nx^2 - 2*Nx +v) = -I;
