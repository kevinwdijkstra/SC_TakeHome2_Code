%% main file
ParameterFile;



%% setup
% % create mesh and matrices
% [D2Mesh, D3Mesh, x, y, z] = CreateMesh(n);
% [D2Mat, D3Mat]   = CreateMatrix(n);
% 
% % create exact solutions
% D2u_ex = u_ex_2D(D2Mesh);
% D3u_ex = u_ex_3D(D3Mesh);
% 
% D2f = f_2D(D2Mesh);
% D3f = f_3D(D3Mesh);

N2D = numel(D2pList);
timeList1 = zeros(N2D,1);
timeList2 = zeros(N2D,1);
NNZ = zeros(N2D,2);
error2D = zeros(N2D,1);

for i = 1:N2D
    disp(strcat("Running p of p max : (",num2str(D2pList(i)),"/",num2str(max(D2pList)),")"))
    %% setup
    n = D2nList(i);
    [D2Mesh, x, y, z] = CreateMesh2D(n);
    [D2Mat]   = CreateMatrix2D(n);
    D2u_ex = u_ex_2D(D2Mesh);

    NNZ(i,1) = nnz(D2Mat);


    %% direct solvers 2D
    % calculate Cholesky Decompostion
    tic;
    C_2D = chol(D2Mat,'lower');
    t_end = toc;
    timeList1(i) = t_end;

    NNZ(i,2) = nnz(C_2D);

    % add boundary conditions
    D2f_dir = CreateBC2D(@(x) u_ex_2D(x),@(x)  f_2D(x),D2Mesh,n);

    
    %% test
    
    % solve 2D problem
    tic;
    u_dir_2D = UpperSolver(C_2D',LowerSolver(C_2D,D2f_dir))';
%     u_dir_2D = C_2D'\(C_2D\(D2f_dir'));
    t_end = toc;
    timeList2(i) = t_end;
    
    error2D(i) = norm(u_dir_2D - D2u_ex',Inf);

end
%% plotting
% loglog(D2nList,mean(timeList1,2))
% hold on
% loglog(D2nList,mean(timeList2,2))
% hold off
% grid on
% legend("Cholesky Decomposition","Solver")


loglog(D2nList,error2D)
grid on
legend("error")



% scatter3(D2Mesh(1,:)',D2Mesh(2,:)',u_dir_2D)
% hold on
% scatter3(D2Mesh(1,:)',D2Mesh(2,:)',D2u_ex)
% hold off
% legend("dir","ex")





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


