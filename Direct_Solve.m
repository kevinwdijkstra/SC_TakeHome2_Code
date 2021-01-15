function [u,times_fac,times_sol,NNZ] = Direct_Solve(A,f,solve_options,n,dim)
% DIRECT_SOLVE Solve the equation Au=f directly using the Cholesky
% decomposition.
%   [u,TIMES_FAC,TIMES_SOL,NNZ] = DIRECT_SOLVE(A,F,SOLVE_OPTIONS,N,DIM)
%   solves AU=F using Cholesky decomposition. It returns a the duration of
%   factorization TIMES_FAC, the duration of solving TIMES_SOL, and the
%   number of nonzero elements NNZ of matrix A and the cholesky factor.

    NNZ = zeros(1,2); % Initialize variable
    NNZ(1) = nnz(A); % Compute the number of nonzero elements of A
    
    % Compute the Cholesky factor of A
    if solve_options.use_MATLAB % Use MATLAB's implementation
        tic;
        C = chol(A,'lower');
        t_end = toc;
    else % Use our implementation
        tic
        C = CholeskyDecompostion(A);
        t_end = toc;
    end
    times_fac = t_end; % Save the duration of factorization
        
    NNZ(2) = nnz(C); % Compute the number of nonzero elements of C

    % Solve the matrix equation
    if solve_options.use_MATLAB % Use MATLAB's implementation
        tic;
        u = C'\(C\(f));
        t_end = toc;
    else % Use our implementation
        tic;
        u = UpperSolver(C',LowerSolver(C,f,(n-1)^(dim-1)),(n-1)^(dim-1));
        t_end = toc;
    end
    times_sol = t_end; % Save the duration of solving
end

