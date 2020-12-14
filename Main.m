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
times_ICBIM_2D = zeros(N2D,1);
times_ICCG_2D = zeros(N2D,1);
NNZ = zeros(N2D,2);
error2D = zeros(N2D,1);
ICBIM_conv_2D = zeros(N2D,M);
ICCG_conv_2D  = zeros(N2D,M);

for i = 1:N2D
    disp(strcat("2D Running p of p max : (",num2str(D2pList(i)),"/",num2str(max(D2pList)),")"))
    %% setup
    n = D2nList(i);
    [D2Mesh, x, y, z] = CreateMesh2D(n);
    [D2Mat]   = CreateMatrix2D(n);
    p = symrcm(D2Mat);
    if ~use_symrcm
        p = 1:length(p);
    end
    D2Mat = D2Mat(p,p);
    D2u_ex = u_ex_2D(D2Mesh);

    
    
    % add boundary conditions
    D2f_dir = CreateBC2D(@(x) u_ex_2D(x),@(x)  f_2D(x),D2Mesh(:,p),n)';


    if do_DIRECT
        %% direct solvers 2D
        fprintf("Direct solver:     ")
        u_dir_2D = 0*D2f_dir;     % we need a blank vector of correct size
        [u_dir_2D(p),times_fac_2D(i),times_sol_2D(i),NNZ(i,:)] = Direct_Solve(D2Mat,D2f_dir,solve_options,n,2);

        error2D(i) = norm(u_dir_2D - D2u_ex',Inf);
        fprintf(strcat(num2str(times_fac_2D(i)+times_sol_2D(i))," seconds\n"));
        
    end
    
    if do_IC_BIM
        %% IC BIM
        fprintf("IC BIM:            ")
        u_k_ICBIM_2D = 0*D2f_dir;     % we need a blank vector of correct size
        [u_k_ICBIM_2D(p),ICBIM_conv_2D(i,:),times_ICBIM_2D(i)] = IC_BIM_Solve(D2Mat,D2f_dir,solve_options,n,2);
        fprintf(strcat(num2str(times_ICBIM_2D(i))," seconds\n"));
    end
    
    if do_ICCG
        %% ICCG
        fprintf("ICCG:              ")
        u_k_ICCG_2D = 0*D2f_dir;     % we need a blank vector of correct size
        [u_k_ICCG_2D(p),ICCG_conv_2D(i,:),times_ICCG_2D(i)] = ICCG_Solve(D2Mat,D2f_dir,solve_options,n,2);
        fprintf(strcat(num2str(times_ICCG_2D(i))," seconds\n"));
    end
    
end


%% 3D problem
N3D = numel(D3pList);
error3D = zeros(N3D,1);
times_fac_3D = zeros(N3D,1);
times_sol_3D = zeros(N3D,1);
times_ICBIM_3D = zeros(N3D,1);
times_ICCG_3D = zeros(N3D,1);
ICBIM_conv_3D = zeros(N3D,M);
ICCG_conv_3D = zeros(N3D,M);

for i = 1:N3D
    disp(strcat("3D Running p of p max : (",num2str(D3pList(i)),"/",num2str(max(D3pList)),")"))
    %% setup
    n = D3nList(i);
    [D3Mesh, x, y, z] = CreateMesh3D(n);
    [D3Mat]   = CreateMatrix3D(n);
    p = symrcm(D3Mat);
    if ~use_symrcm
        p = 1:length(p);
    end
    D3Mat = D3Mat(p,p);
    D3u_ex = u_ex_3D(D3Mesh);
    
    % add boundary conditions
    D3f_dir = CreateBC3D(@(x) u_ex_3D(x),@(x)  f_3D(x),D3Mesh(:,p),n)';

    if do_DIRECT
        %% direct solvers 3D
        fprintf("Direct solver:     ")

        u_dir_3D = 0*D3f_dir;     % we need a blank vector of correct size
        [u_dir_3D(p),times_fac_3D(i),times_sol_3D(i),~] = Direct_Solve(D3Mat,D3f_dir,solve_options,n,3);

        error3D(i) = norm(u_dir_3D - D3u_ex',Inf);
        fprintf(strcat(num2str(times_fac_3D(i)+times_sol_3D(i))," seconds\n"));
    end
    
    if do_IC_BIM
        %% IC BIM
        fprintf("IC BIM:            ")
        u_k_ICBIM_3D = 0*D3f_dir;     % we need a blank vector of correct size
        [u_k_ICBIM_3D(p),ICBIM_conv_3D(i,:),times_ICBIM_3D(i)] = IC_BIM_Solve(D3Mat,D3f_dir,solve_options,n,3);
        fprintf(strcat(num2str(times_ICBIM_3D(i))," seconds\n"));
    end
    
    if do_ICCG
        %% ICCG
        fprintf("ICCG:              ")
        u_k_ICCG_3D = 0*D3f_dir;     % we need a blank vector of correct size
        [u_k_ICCG_3D(p),ICCG_conv_3D(i,:),times_ICCG_3D(i)] = ICCG_Solve(D3Mat,D3f_dir,solve_options,n,3);
        fprintf(strcat(num2str(times_ICCG_3D(i))," seconds\n"));
    end
    
end

%% plotting

if plot_figure
    disp("Plotting")
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

    subplot(4,1,1)
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
    title("2D Direct:");

    subplot(4,1,2)
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
    title("3D Direct:");

    subplot(4,1,3)
    loglog(D2nList,times_ICBIM_2D);
    hold on
    loglog(D3nList,times_ICBIM_3D);
    hold off
    set(gca,'xtick',D2nList);
    set (gca, 'XTickLabel', strcat('2^{',num2str((D2pList(:))),'}'));
    grid on
    legend("2D","3D",'location','southeast');
    xlabel("number of grid elements in each dimension");
    ylabel("time in seconds of each operation");
    title("IC BIM:");
    
    subplot(4,1,4)
    loglog(D2nList,times_ICCG_2D);
    hold on
    loglog(D3nList,times_ICCG_3D);
    hold off
    set(gca,'xtick',D2nList);
    set (gca, 'XTickLabel', strcat('2^{',num2str((D2pList(:))),'}'));
    grid on
    legend("2D","3D",'location','southeast');
    xlabel("number of grid elements in each dimension");
    ylabel("time in seconds of each operation");
    title("ICCG:");


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
    legend(compose('p=%u',D2pList))
    subplot(2,1,2)
    loglog(ICBIM_conv_3D');
    title("3D IC BIM convergence")
    grid on
    xlabel("Number of iterations");
    ylabel("2 norm of the error");
    legend(compose('p=%u',D3pList))

    %% convergence of ICCG
    figure(5)
    subplot(2,1,1)
    loglog(ICCG_conv_2D');
    title("2D ICCG convergence")
    grid on
    xlabel("Number of iterations");
    ylabel("2 norm of the error");
    legend(compose('p=%u',D2pList))
    subplot(2,1,2)
    loglog(ICCG_conv_3D');
    title("3D ICCG convergence")
    grid on
    xlabel("Number of iterations");
    ylabel("2 norm of the error");
    legend(compose('p=%u',D3pList))
    
    %% total excecution time
    figure(6)
    subplot(2,1,1)
    loglog(D2nList,times_sol_2D+times_fac_2D);
    hold on
    loglog(D2nList,times_ICBIM_2D);
    loglog(D2nList,times_ICCG_2D);
    hold off
    title("time of convergence in 2D")
    grid on
    xlabel("number of grid elements in each dimension");
    ylabel("total time");
    legend("Direct","IC BIM","ICCG",'location','northwest')
    subplot(2,1,2)
    loglog(D3nList,times_sol_3D+times_fac_3D);
    hold on
    loglog(D3nList,times_ICBIM_3D);
    loglog(D3nList,times_ICCG_3D);
    hold off
    title("time of convergence in 3D")
    grid on
    xlabel("number of grid elements in each dimension");
    ylabel("total time");
    legend("Direct","IC BIM","ICCG",'location','northwest')
    
    
    set(0,'DefaultFigureWindowStyle','normal')
end

disp("Done");





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


