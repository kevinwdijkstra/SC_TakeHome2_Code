%% main file
p = 5;
N = 2^p;



%% setup
% create mesh and matrices
[D2Mesh, D3Mesh, x, y, z] = CreateMesh(N);
[D2Mat, D3Mat]   = CreateMatrix(N);

% create exact solutions
D2u_ex = u_ex_2D(D2Mesh);
D3u_ex = u_ex_3D(D3Mesh);

D2f = f_2D(D2Mesh);
D3f = f_3D(D3Mesh);


%% direct solvers
% 2D
[L_2D,U_2D] = lu(D2Mat);

% add boundary conditions
Ixlow = (D2Mesh(1,:)-1/(N)==0);
Ixhig = (D2Mesh(1,:)+1/(N)==1);
Iylow = (D2Mesh(2,:)-1/(N)==0);
Iyhig = (D2Mesh(2,:)+1/(N)==1);

D2f_dir     = D2f + n^2*(   Ixlow.*u_ex_2D(D2Mesh - [1/(n);0])+...
                            Ixhig.*u_ex_2D(D2Mesh + [1/(n);0])+...
                            Iylow.*u_ex_2D(D2Mesh - [0;1/(n)])+...
                            Iyhig.*u_ex_2D(D2Mesh + [0;1/(n)]));

u_dir_2D = D2Mat\D2f_dir';


scatter3(D2Mesh(1,:)',D2Mesh(2,:)',u_dir_2D)
hold on
scatter3(D2Mesh(1,:)',D2Mesh(2,:)',D2u_ex)
hold off
legend("dir","ex")













%% equations
function u = u_ex_2D(x)
    u = x(1,:).^4.*x(2,:).^5;
end

function u = u_ex_3D(x)
    u = x(1,:).^4.*x(2,:).^5.*x(3,:).^6;
end

function u = f_2D(x)
    u = -12*x(1,:).^2.*x(2,:).^5-20*x(1,:).^4.*x(2,:).^3;
end

function u = f_3D(x)
    u = -12*x(1,:).^2.*x(2,:).^5.*x(3,:).^6-20*x(1,:).^4.*x(2,:).^3.*x(3,:).^6-30*x(1,:).^4.*x(2,:).^5.*x(3,:).^4;
end


