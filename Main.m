%% main file
ParameterFile;



%% setup
% create mesh and matrices
[D2Mesh, D3Mesh, x, y, z] = CreateMesh(n);
[D2Mat, D3Mat]   = CreateMatrix(n);

% create exact solutions
D2u_ex = u_ex_2D(D2Mesh);
D3u_ex = u_ex_3D(D3Mesh);

D2f = f_2D(D2Mesh);
D3f = f_3D(D3Mesh);

N2D = numel(D2pList);
timeList = zeros(N,2);

for i = 1:N2D
    %% direct solvers 2D
    % calculate Cholesky Decompostion
    C_2D = chol(D2Mat);

    % add boundary conditions
    D2f_dir = CreateBC2D(@(x) u_ex_2D(x),@(x)  f_2D(x),D2Mesh,n);

    % solve 2D problem
    u_dir_2D = C_2D\(C_2D'\D2f_dir');

end




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


