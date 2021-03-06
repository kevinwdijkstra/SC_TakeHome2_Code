%% main file
ParameterFile;



%% setup

N2D = numel(D2pList);               % number of different 2D simulations to run

% pre alocate memory
times_fac_2D = zeros(N2D,1);        
times_sol_2D = zeros(N2D,1);
times_ICBIM_2D = zeros(N2D,1);
times_ICCG_2D = zeros(N2D,1);
NNZ = zeros(N2D,2);
error2D = zeros(N2D,1);
ICBIM_conv_2D = zeros(N2D,M); % M is an upper bound of total iterations
ICCG_conv_2D  = zeros(N2D,M); % M is an upper bound of total iterations

for i = 1:N2D % loop over different 2D simulations
    disp(strcat("2D Running p of p max : (",num2str(D2pList(i)),"/",num2str(max(D2pList)),")"))
    %% setup of mesh/matrix and right hand side vector.
    n = D2nList(i);
    [D2Mesh, x, y, z] = CreateMesh2D(n);
    [D2Mat]   = CreateMatrix2D(n);
    if use_symrcm % wether to use the reordering scheme for the matrix
        p = symrcm(D2Mat); % p is an index vector that describes the reordering
    else
        p = 1:(n-1)^2;
    end
    % calculate exact solution
    D2u_ex = u_ex_2D(D2Mesh);

    % setup right hand side vector
    D2f_dir = CreateBC2D(@(x) u_ex_2D(x),@(x)  f_2D(x),D2Mesh,n)';


    if do_DIRECT
        %% direct solvers 2D
        fprintf("Direct solver:     ")
        % we need a blank vector of correct size to store result
        u_dir_2D = 0*D2f_dir;     
        % run direct solver
        [u_dir_2D(p),times_fac_2D(i),times_sol_2D(i),NNZ(i,:)] = Direct_Solve(D2Mat(p,p),D2f_dir(p),solve_options,n,2); 
        
        % calculate error
        error2D(i) = norm(u_dir_2D - D2u_ex',Inf);
        % print total execution time
        fprintf(strcat(num2str(times_fac_2D(i)+times_sol_2D(i))," seconds\n"));
        
    end
    
    if do_IC_BIM
        %% IC BIM
        fprintf("IC BIM:            ")
        % we need a blank vector of correct size
        u_k_ICBIM_2D = 0*D2f_dir;     
        % run IC BIM solver
        [u_k_ICBIM_2D,ICBIM_conv_2D(i,:),times_ICBIM_2D(i)] = IC_BIM_Solve(D2Mat,D2f_dir,solve_options,n,2); 
        % print total execution time
        fprintf(strcat(num2str(times_ICBIM_2D(i))," seconds\n"));
    end
    
    if do_ICCG
        %% ICCG
        fprintf("ICCG:              ")
        % we need a blank vector of correct size
        u_k_ICCG_2D = 0*D2f_dir;  
        % run ICCG solver
        [u_k_ICCG_2D,ICCG_conv_2D(i,:),times_ICCG_2D(i)] = ICCG_Solve(D2Mat,D2f_dir,solve_options,n,2);
        % print total execution time
        fprintf(strcat(num2str(times_ICCG_2D(i))," seconds\n"));
    end
    
end


%% 3D problem
% set total number of 3D simulations to run
N3D = numel(D3pList);
% pre allocate memory
error3D = zeros(N3D,1);
times_fac_3D = zeros(N3D,1);
times_sol_3D = zeros(N3D,1);
times_ICBIM_3D = zeros(N3D,1);
times_ICCG_3D = zeros(N3D,1);
ICBIM_conv_3D = zeros(N3D,M); % M is an upper bound of total iterations
ICCG_conv_3D = zeros(N3D,M);  % M is an upper bound of total iterations

for i = 1:N3D
    disp(strcat("3D Running p of p max : (",num2str(D3pList(i)),"/",num2str(max(D3pList)),")"))
    %%  setup of mesh/matrix and right hand side vector.
    n = D3nList(i);
    [D3Mesh, x, y, z] = CreateMesh3D(n);
    [D3Mat]   = CreateMatrix3D(n);
    if use_symrcm % wether te use the reordering scheme
        p = symrcm(D3Mat); % p is an index vector that describes the reordering
    else
        p = 1:(n-1)^3;
    end
    % calculate exact solution
    D3u_ex = u_ex_3D(D3Mesh);
    
    % create right hand side vector
    D3f_dir = CreateBC3D(@(x) u_ex_3D(x),@(x)  f_3D(x),D3Mesh,n)';

    if do_DIRECT
        %% direct solvers 3D
        fprintf("Direct solver:     ")
        % we need a blank vector of correct size
        u_dir_3D = 0*D3f_dir;     
        % do direct solver
        [u_dir_3D(p),times_fac_3D(i),times_sol_3D(i),~] = Direct_Solve(D3Mat(p,p),D3f_dir(p),solve_options,n,3);

        % calculate error
        error3D(i) = norm(u_dir_3D - D3u_ex',Inf);
        % print total execution time
        fprintf(strcat(num2str(times_fac_3D(i)+times_sol_3D(i))," seconds\n"));
    end
    
    if do_IC_BIM
        %% IC BIM
        fprintf("IC BIM:            ")
        % we need a blank vector of correct size
        u_k_ICBIM_3D = 0*D3f_dir;     
        % do IC BIM solve
        [u_k_ICBIM_3D,ICBIM_conv_3D(i,:),times_ICBIM_3D(i)] = IC_BIM_Solve(D3Mat,D3f_dir,solve_options,n,3);
        % print total execution time
        fprintf(strcat(num2str(times_ICBIM_3D(i))," seconds\n"));
    end
    
    if do_ICCG
        %% ICCG
        fprintf("ICCG:              ")
        % we need a blank vector of correct size
        u_k_ICCG_3D = 0*D3f_dir;     
        % do ICCG solve
        [u_k_ICCG_3D,ICCG_conv_3D(i,:),times_ICCG_3D(i)] = ICCG_Solve(D3Mat,D3f_dir,solve_options,n,3);
        % print total execution time
        fprintf(strcat(num2str(times_ICCG_3D(i))," seconds\n"));
    end
    
end

%% plotting

if plot_figure
    disp("Plotting")
    set(0,'DefaultFigureWindowStyle','docked')

  
    %% error solutions
    figure(1)
    
    
%     b=axes('Position',[.1 .1 .8 1e-12]);
    
    
    % axis for m/s
    b=axes('Position',[.1 .1 .8 1e-15]);
    set(b,'Units','normalized');
    set(b,'Color','none');
    
    loglog((D2nList),error2D,'k','LineWidth',0.01)
    set(b,'xtick',D2nList);
    set(b, 'XTickLabel', strcat('2^{',num2str((D2pList(:))),'}'));

    a=axes('Position',[.1 .1 .8 .8]);
    loglog((D2nList),error2D)
    hold on
    loglog((D3nList),error3D)
    hold off
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
    ylabel("time [seconds]");
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
    ylabel("time [seconds]");
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
    ylabel("time [seconds]");
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
    ylabel("time [seconds]");
    title("ICCG:");


    %% NNZ fill in ratio
    figure(3)
    plot(D2nList,NNZ(:,2)./NNZ(:,1),'--o');
%     set(gca,'xtick',D2nList);
%     set (gca, 'XTickLabel', strcat('2^{',num2str((D2pList(:))),'}'));
    grid on
    xlabel("number of grid elements in each dimension");
    ylabel("$\frac{NNZ(C_h)}{NNZ(A_h)}$",'interpreter','latex');
    title("the $\frac{NNZ(C_h)}{NNZ(A_h)}$ ratio for different matrix sizes.",'interpreter','latex');

    %% convergence of IC BIM
    figure(4)
    subplot(2,1,1)
    loglog(ICBIM_conv_2D');
    ylim([1e-10 1])
    title("2D IC BIM convergence")
    grid on
    xlabel("Number of iterations");
    ylabel("2 norm of the error");
    legend(compose('p=%u',D2pList))
    subplot(2,1,2)
    loglog(ICBIM_conv_3D');
    ylim([1e-10 1])
    title("3D IC BIM convergence")
    grid on
    xlabel("Number of iterations");
    ylabel("2 norm of the error");
    legend(compose('p=%u',D3pList))

    %% convergence of ICCG
    figure(5)
    subplot(2,1,1)
    loglog(ICCG_conv_2D');
    ylim([1e-10 1])
    title("2D ICCG convergence")
    grid on
    xlabel("Number of iterations");
    ylabel("2 norm of the error");
    legend(compose('p=%u',D2pList))
    subplot(2,1,2)
    loglog(ICCG_conv_3D');
    ylim([1e-10 1])
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
    
    
    %% Last 5 iterations convergence ratio
    if do_IC_BIM
    figure(7)
    % 2D
    subplot(2,1,1)
    for i = 1:numel(D2pList)
        I2D = find(ICBIM_conv_2D(i,:) > 0);
        Iend = I2D(end);
        plot([-5 -4 -3 -2 -1 0],ICBIM_conv_2D(i,Iend-5:Iend)./(ICBIM_conv_2D(i,Iend-6:Iend-1)),"-o");
        hold on
    end
    hold off
    legend(num2str(D2pList'))
    grid on
    title("2D last 5 iterations ratio of convergence")
    xlabel("iteration n from last iteration")
    ylabel("ratio $\frac{||\mathbf{r}_n||_2}{||\mathbf{r}_{n-1}||_2}$",'interpreter','latex')
    
    % 3D
    subplot(2,1,2)
    for i = 1:numel(D3pList)
        I3D = find(ICBIM_conv_3D(i,:) > 0);
        Iend = I3D(end);
        plot([-5 -4 -3 -2 -1 0],ICBIM_conv_3D(i,Iend-5:Iend)./(ICBIM_conv_3D(i,Iend-6:Iend-1)),"-o");
        hold on
    end
    hold off
    legend(num2str(D3pList'))
    grid on
    title("3D last 5 iterations ratio of convergence")
    xlabel("iteration n from last iteration")
    ylabel("ratio $\frac{||\mathbf{r}_n||_2}{||\mathbf{r}_{n-1}||_2}$",'interpreter','latex')
    end
    
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


