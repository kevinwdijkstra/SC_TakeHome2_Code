function [uk,ICBIM_conv,times_IC] = IC_BIM_Solve(A,f,solve_options,n,dim)
% IC_BIM_SOLVE Solve the equation Au=f directly using the Incomplete
% Cholesky decomposition as a Basic Iterative Method.
%   [UK,ICBIM_CONV,TIMES_IC] = IC_BIM_SOLVE(A,F,SOLVE_OPTIONS,N,DIM)
%   solves AU=F using the Incomplete Cholesky decomposition as BIM. It
%   returns the solution UK, the required time before convergence TIMES_IC
%   and the residuals of each step ICBIM_CONV.

    ICBIM_conv = zeros(1,solve_options.M); % Initialize vector of residuals

    uk = 0*f; % initialize solution vector
    rk = f; % initialize residuals vector
    crit = solve_options.epsilon*norm(f); % set stopping criterion for rk
    j = 0;
    nF = norm(f); % save norm of f
    tic;
    if solve_options.use_MATLAB % use MATLAB's implementations
        % Here the splitting is used A = L*L' - R
        L = ichol(A); % Compute IC matrix
        R = L*L' - A;
        while norm(rk)>crit % While stopping criterion not reached
            uk = L'\(L\(R*uk + f)); % Compute next solution
            rk = f - A*uk; % Compute residual
            j = j+1; % increment
            ICBIM_conv(1,j) = norm(rk)/nF; % Save scaled residual
        end
    else % Use our implementations
        % Here the splitting is used A = L*D^-1*L' - R
        L = IncompleteCholesky(A,n-1,dim); % Compute IC matrix
        Dinv = inv(spdiags(spdiags(L,0),0,(n-1)^dim,(n-1)^dim)); % Compute D^-1
        L1 = L*Dinv; % Precompute lower triangular matrix
        R = L1*L' - A; % Compute matrix R
        while norm(rk)>crit % While stopping criterion not reached
            uk = UpperSolver(L',LowerSolver(L1,R*uk + f,(n-1)^(dim-1)),...
                (n-1)^(dim-1)); % Compute next solution
            rk = f - A*uk; % Compute residuals
            j = j+1; % increment
            ICBIM_conv(1,j) = norm(rk)/nF; % Save scaled residual
        end
    end
    times_IC = toc;
end

