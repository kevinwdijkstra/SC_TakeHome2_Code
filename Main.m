%% main file
ParameterFile;

%% still to do
% 3D



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
times_fac_2D = zeros(N2D,1);
times_sol_2D = zeros(N2D,1);
times_IC_2D = zeros(N2D,1);
NNZ = zeros(N2D,2);
error2D = zeros(N2D,1);
ICBIM_conv_2D = zeros(N2D,M);

for i = 1:N2D
    disp(strcat("2D Running p of p max : (",num2str(D2pList(i)),"/",num2str(max(D2pList)),")"))
    %% setup
    n = D2nList(i);
    [D2Mesh, x, y, z] = CreateMesh2D(n);
    [D2Mat]   = CreateMatrix2D(n);
    D2u_ex = u_ex_2D(D2Mesh);

    NNZ(i,1) = nnz(D2Mat);
    
    % add boundary conditions
    D2f_dir = CreateBC2D(@(x) u_ex_2D(x),@(x)  f_2D(x),D2Mesh,n)';


    %% direct solvers 2D
    disp("Direct solver")
    % calculate Cholesky Decompostion
    if use_MATLAB
        tic;
        C_2D = chol(D2Mat,'lower');
        t_end = toc;
    else
        tic
        C_2D = CholeskyDecompostion(D2Mat);
        t_end = toc;
    end
    times_fac_2D(i) = t_end;
        
%     C_2D1 = chol(D2Mat,'lower');
%     C_2D2 = CholeskyDecompostion(D2Mat);
%     error = sum(sum(abs(C_2D1-C_2D2)))
   
    NNZ(i,2) = nnz(C_2D);

    % solve 2D problem
    if use_MATLAB
        tic;
        u_dir_2D = C_2D'\(C_2D\(D2f_dir));
        t_end = toc;
    else
        tic;
        u_dir_2D = UpperSolver(C_2D',LowerSolver(C_2D,D2f_dir))';
        t_end = toc;
    end
    times_sol_2D(i) = t_end;
    
    error2D(i) = norm(u_dir_2D - D2u_ex',Inf);
    
    
    %% IC BIM
    disp("IC BIM")
    uk = 0*D2f_dir;
    rk = D2f_dir;
    crit = epsilon*norm(D2f_dir);
    j = 0;
    tic;
    if use_MATLAB
        L = ichol(D2Mat);
        R = L*L' - D2Mat;
        while norm(rk)>crit
            uk = L'\(L\(R*uk + D2f_dir));
            rk = D2f_dir - D2Mat*uk;
            j = j+1;
            ICBIM_conv_2D(i,j) = norm(rk)/norm(D2f_dir);
        end
    else
        L = IncompleteCholesky(D2Mat,n-1);
        Dinv = inv(spdiags(spdiags(L,0),0,(n-1)^2,(n-1)^2));
        R = L*Dinv*L' - D2Mat;
        while norm(rk)>crit
            uk = UpperSolver(L',LowerSolver(L*Dinv,R*uk + D2f_dir))';
            rk = D2f_dir - D2Mat*uk;
            j = j+1;
            ICBIM_conv_2D(i,j) = norm(rk)/norm(D2f_dir);
        end
    end
    times_IC_2D(i) = toc;
    

end


%% 3D problem
N3D = numel(D3pList);
error3D = zeros(N3D,1);
times_fac_3D = zeros(N3D,1);
times_sol_3D = zeros(N3D,1);
times_IC_3D = zeros(N3D,1);
ICBIM_conv_3D = zeros(N3D,M);

for i = 1:N3D
    disp(strcat("3D Running p of p max : (",num2str(D3pList(i)),"/",num2str(max(D3pList)),")"))
    %% setup
    n = D3nList(i);
    [D3Mesh, x, y, z] = CreateMesh3D(n);
    [D3Mat]   = CreateMatrix3D(n);
    D3u_ex = u_ex_3D(D3Mesh);
    
    % add boundary conditions
    D3f_dir = CreateBC3D(@(x) u_ex_3D(x),@(x)  f_3D(x),D3Mesh,n)';


    %% direct solvers 3D
    % calculate Cholesky Decompostion
    disp("Direct Solver")

    if use_MATLAB
        tic;
        C_3D = chol(D3Mat,'lower');
        times_fac_3D(i) = toc;
    else
        tic;
        C_3D = CholeskyDecompostion(D3Mat);
        times_fac_3D(i) = toc;
    end
        
    % solve 3D problem
    if use_MATLAB
        tic;
        u_dir_3D = C_3D'\(C_3D\(D3f_dir));
        times_sol_3D(i) = toc;
    else
        tic;
        u_dir_3D = UpperSolver(C_3D',LowerSolver(C_3D,D3f_dir))';
        times_sol_3D(i) = toc;
    end
    
    error3D(i) = norm(u_dir_3D - D3u_ex',Inf);
    
        %% IC BIM
    disp("IC BIM")
    uk = 0*D3f_dir;
    rk = D3f_dir;
    crit = epsilon*norm(D3f_dir);
    j = 0;
    tic;
    if use_MATLAB
        L = ichol(D3Mat);
        R = L*L' - D3Mat;
        while norm(rk)>crit
            uk = L'\(L\(R*uk + D3f_dir));
            rk = D3f_dir - D3Mat*uk;
            j = j+1;
            ICBIM_conv_3D(i,j) = norm(rk)/norm(D3f_dir);
        end
    else
        L = IncompleteCholesky(D3Mat,[n-1,(n-1)^2]);
        Dinv = inv(spdiags(spdiags(L,0),0,(n-1)^3,(n-1)^3));
        L1 = L*Dinv;
        R = L*Dinv*L' - D3Mat;
        while norm(rk)>crit
            uk = UpperSolver(L',LowerSolver(L1,R*uk + D3f_dir))';
            rk = D3f_dir - D3Mat*uk;
            j = j+1;
            ICBIM_conv_3D(i,j) = norm(rk)/norm(D3f_dir);
        end
    end
    times_IC_3D(i) = toc;

end

%% plotting
set(0,'DefaultFigureWindowStyle','docked')

%% error solutions
figure(1)
loglog((D2nList),error2D)
hold on
loglog((D3nList),error3D)
hold off
set(gca,'xtick',D2nList);
set (gca, 'XTickLabel', strcat('2^{',num2str((D2pList(:))),'}'));
grid on
legend("2D","3D")
xlabel("number of grid elements in each dimension")
ylabel("inf norm of the error")
title(["Error convergence of the discritized system in 2D and 3D";" with the help of direct solvers for various number of grid spacings."]);

%% time for operations
figure(2)
title({'Time requirements of the Cholesky decomposition and';'the Forward/Backward solve steps for various number of grid elements.'});

subplot(3,1,1)
loglog(D2nList,times_fac_2D);
hold on
loglog(D2nList,times_sol_2D);
hold off
set(gca,'xtick',D2nList);
set (gca, 'XTickLabel', strcat('2^{',num2str((D2pList(:))),'}'));
grid on
legend("Cholesky decomposition","Forward/Backward solve",'location','southeast');
xlabel("number of grid elements in each dimension");
ylabel("time in seconds of each operation");

subplot(3,1,2)
loglog(D3nList,times_fac_3D);
hold on
loglog(D3nList,times_sol_3D);
hold off
set(gca,'xtick',D3nList);
set (gca, 'XTickLabel', strcat('2^{',num2str((D3pList(:))),'}'));
grid on
legend("Cholesky decomposition","Forward/Backward solve",'location','southeast');
xlabel("number of grid elements in each dimension");
ylabel("time in seconds of each operation");

subplot(3,1,3)
loglog(D2nList,times_IC_2D);
hold on
loglog(D3nList,times_IC_3D);
hold off
set(gca,'xtick',D2nList);
set (gca, 'XTickLabel', strcat('2^{',num2str((D2pList(:))),'}'));
grid on
legend("Cholesky decomposition","Forward/Backward solve",'location','southeast');
xlabel("number of grid elements in each dimension");
ylabel("time in seconds of each operation");


%% NNZ fill in ratio
figure(3)
semilogx(D2nList,NNZ(:,2)./NNZ(:,1));
set(gca,'xtick',D2nList);
set (gca, 'XTickLabel', strcat('2^{',num2str((D2pList(:))),'}'));
grid on
xlabel("number of grid elements in each dimension");
ylabel("ratio of NNZ of C_h and A_h");
title("the C_h/A_h ratio of the non zero elements counts in each matrix.");

%% convergence of IC BIM
figure(4)
subplot(2,1,1)
loglog(ICBIM_conv_2D');
title("2D IC BIM convergence")
grid on
xlabel("Number of iterations");
ylabel("2 norm of the error");
legend(compose('p=%u',D3pList))
subplot(2,1,2)
loglog(ICBIM_conv_3D');
title("3D IC BIM convergence")
grid on
xlabel("Number of iterations");
ylabel("2 norm of the error");
legend(compose('p=%u',D3pList))





set(0,'DefaultFigureWindowStyle','normal')



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
    u = -12*x(1,:).^2.*x(2,:).^5.*x(3,:).^6 + ...
        -20*x(1,:).^4.*x(2,:).^3.*x(3,:).^6 + ...
        -30*x(1,:).^4.*x(2,:).^5.*x(3,:).^4;
end


