function [u_k,ICCG_conv,times_ICCG] = ICCG_Solve(A,f,solve_options,n,dim)
%% solves Au = f with the help of ICCG
% inputs:
% - A problem matrix
% - f right hand side vector
% - solve_options contains various settings, e.g. stopping criterion epsilon
% - n number of elements in each dimension
% - dim number of dimensions

% preallocate memory in struct convergence to be passed between functions 
convergence.ICCG_conv = zeros(1,solve_options.M);
convergence.i = 1;
convergence.nF = norm(f);
crit = solve_options.epsilon*norm(f);

% calculate preconditioner
if solve_options.use_MATLAB % wether to use matlab or our implementation
    M_pre = ichol(A);
else
    M_pre = IncompleteCholesky(A,n-1,dim);
end    

% set initial values to zero
u_k = 0*f;        
r_k = f;

% for timing purposes
tic 


if solve_options.use_MATLAB % wether to use matlab or our implementation
    % do k = 1 as the code is a little different
    % calculate pre conditioning
    z_k = M_pre'\(M_pre\r_k);
    
    p_k = z_k;                          % set p_k
    A_p_k = A*p_k;                  	% pre calculate Ap_{k}
    n_r_k_old = r_k'*z_k;               % calculate r_{k}*z_{k}
    alpha_k = n_r_k_old/(p_k'*A_p_k);   % calculate alpha with the help of previously calculated values
    u_k = u_k + alpha_k*p_k;            % update iteratate
    r_k = r_k - alpha_k*A_p_k;          % update residual

    % r_k might differ from this criterion due to rounding errors. We start
    % again if that is the case
    while norm(f-A*u_k)>crit
        % do loop
        [u_k,convergence] = while_loop_MATLAB(u_k,r_k,p_k,n_r_k_old,crit,A,M_pre,convergence);
    end
else
    % do k = 1 as the code is a little different
    % calculate pre conditioning
    z_k = UpperSolver(M_pre',LowerSolver(M_pre,r_k,(n-1)^(dim-1)),(n-1)^(dim-1));
    
    p_k = z_k;                          % set p_k
    A_p_k = A*p_k;                      % pre calculate Ap_{k}
    n_r_k_old = r_k'*z_k;               % calculate r_{k}*z_{k}
    alpha_k = n_r_k_old/(p_k'*A_p_k);   % calculate alpha with the help of previously calculated values
    u_k = u_k + alpha_k*p_k;            % update iteratate
    r_k = r_k - alpha_k*A_p_k;          % update residual

    % r_k might differ from this criterion due to rounding errors. We start
    % again if that is the case
    while norm(f-A*u_k)>crit
        [u_k,convergence] = while_loop(u_k,r_k,p_k,n_r_k_old,crit,A,M_pre,convergence,n,dim);
    end
end  

% store total execution time
times_ICCG = toc;

% get iteration values from convergence struct.
ICCG_conv = convergence.ICCG_conv;
    
end

function [u_k,convergence] = while_loop_MATLAB(u_k,r_k,p_k,n_r_k_old,crit,A,M_pre,convergence)
    % this functions uses the solvers of MATLAB
    
    % loop till the norm satisfies the criterium, note that this is
    % epsilon*norm(f).
    while norm(r_k)>crit
        % calculate pre conditioning
        z_k = M_pre'\(M_pre\r_k);
        
        n_r_k = r_k'*z_k;               % calculate r_{k-1}*z_{k-1}
        beta_k = n_r_k/n_r_k_old;       % calculate (r_{k-1}*z_{k-1})/(r_{k-2}*z_{k-2}) by using previously calculated values
        p_k = z_k + beta_k*p_k;         % calculate p_{k} from previous values
        
        A_p_k = A*p_k;                  % pre calculate Ap_{k}
        alpha_k = n_r_k/(p_k'*A_p_k);   % calculate alpha with the help of previously calculated values
        u_k = u_k + alpha_k*p_k;        % update iteratate
        r_k = r_k - alpha_k*A_p_k;      % update residual
        
        % store current iteration for plotting
        convergence.ICCG_conv(convergence.i) = norm(r_k)/convergence.nF;
        convergence.i = convergence.i + 1;
        
        % updatate old inner product with new value for next loop.
        n_r_k_old = n_r_k;              % k-2 <= k-1
    end
end

function [u_k,convergence] = while_loop(u_k,r_k,p_k,n_r_k_old,crit,A,M_pre,convergence,n,dim)
    % this functions uses the solvers we wrote
    
    % loop till the norm satisfies the criterium, note that this is
    % epsilon*norm(f).
    while norm(r_k)>crit
        % calculate pre conditioning
        z_k = UpperSolver(M_pre',LowerSolver(M_pre,r_k,(n-1)^(dim-1)),(n-1)^(dim-1));
        
        n_r_k = r_k'*z_k;               % calculate r_{k-1}*z_{k-1}
        beta_k = n_r_k/n_r_k_old;       % calculate (r_{k-1}*z_{k-1})/(r_{k-2}*z_{k-2}) by using previously calculated values
        p_k = z_k + beta_k*p_k;         % calculate p_{k} from previous values
        
        A_p_k = A*p_k;              	% pre calculate Ap_{k}
        alpha_k = n_r_k/(p_k'*A_p_k);   % calculate alpha with the help of previously calculated values
        u_k = u_k + alpha_k*p_k;        % update iteratate
        r_k = r_k - alpha_k*A_p_k;      % update residual
        
        % store current iteration for plotting
        convergence.ICCG_conv(convergence.i) = norm(r_k)/convergence.nF;
        convergence.i = convergence.i + 1;
        
        % updatate old inner product with new value for next loop.
        n_r_k_old = n_r_k;              % k-2 <= k-1
    end
end

