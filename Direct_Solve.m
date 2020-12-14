function [u,times_fac,times_sol,NNZ] = Direct_Solve(A,f,solve_options)

    NNZ = zeros(1,2);
    NNZ(1) = nnz(A);
    if solve_options.use_MATLAB
        tic;
        C = chol(A,'lower');
        t_end = toc;
    else
        tic
        C = CholeskyDecompostion(A);
        t_end = toc;
    end
    times_fac = t_end;
        
   
    NNZ(2) = nnz(C);

    % solve 2D problem
    if solve_options.use_MATLAB
        tic;
        u = C'\(C\(f));
        t_end = toc;
    else
        tic;
        u = UpperSolver(C',LowerSolver(C,f,size(A,1)-1),size(A,1)-1);
        t_end = toc;
    end
    times_sol = t_end;
    % reverse reordering
    

end

