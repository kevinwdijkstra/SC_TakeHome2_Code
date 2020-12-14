function [u_k,ICCG_conv,times_ICCG] = ICCG_Solve(A,f,solve_options)

    convergence.ICCG_conv = zeros(1,solve_options.M);
    convergence.i = 1;
    convergence.nF = norm(f);
    crit = solve_options.epsilon*norm(f);
    
    % calculate preconditioner
    if solve_options.use_MATLAB
        M_pre = ichol(A);
    else
        M_pre = IncompleteCholesky(A,size(A,1)-1);
    end    
    
    % k =0
    u_k = 0*f;        
    r_k = f;
   
    tic 
    
    % do k = 1
    if solve_options.use_MATLAB
        z_k = M_pre'\(M_pre\r_k);
        p_k = z_k;
        A_p_k = A*p_k;                  % k
        n_r_k_old = r_k'*z_k;
        alpha_k = n_r_k_old/(p_k'*A_p_k);   % k-1 / k
        u_k = u_k + alpha_k*p_k;            % k
        r_k = r_k - alpha_k*A_p_k;          % k - k-1
    else
        z_k = UpperSolver(M_pre',LowerSolver(M_pre,r_k,size(A,1)-1),size(A,1)-1);
        p_k = z_k;
        A_p_k = A*p_k;                  % k
        n_r_k_old = r_k'*z_k;
        alpha_k = n_r_k_old/(p_k'*A_p_k);   % k-1 / k
        u_k = u_k + alpha_k*p_k;            % k
        r_k = r_k - alpha_k*A_p_k;          % k - k-1
    end  


    
    % r_k might differ from this criterion. We start again if that is the
    % case
    while norm(f-A*u_k)>crit
        [u_k,convergence] = while_loop(u_k,r_k,p_k,n_r_k_old,crit,A,M_pre,convergence);
    end
    times_ICCG = toc;
    
    ICCG_conv = convergence.ICCG_conv;
    
    
end

function [u_k,convergence] = while_loop(u_k,r_k,p_k,n_r_k_old,crit,A,M_pre,convergence)
    while norm(r_k)>crit
        z_k = M_pre'\(M_pre\r_k);
        n_r_k = r_k'*z_k;               % k-1
        beta_k = n_r_k/n_r_k_old;       % k-1 / k-2
        p_k = z_k + beta_k*p_k;         % k-1
        
        A_p_k = A*p_k;              % k
        alpha_k = n_r_k/(p_k'*A_p_k);   % k-1 / k
        u_k = u_k + alpha_k*p_k;        % k
        r_k = r_k - alpha_k*A_p_k;      % k - k-1
        convergence.ICCG_conv(convergence.i) = norm(r_k)/convergence.nF;
        convergence.i = convergence.i + 1;
        
        n_r_k_old = n_r_k;              % k-2 <= k-1
    end
end

