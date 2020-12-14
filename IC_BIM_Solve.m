function [uk,ICBIM_conv,times_IC] = IC_BIM_Solve(A,f,solve_options,n,dim)


    ICBIM_conv = zeros(1,solve_options.M);

    uk = 0*f;
    rk = f;
    crit = solve_options.epsilon*norm(f);
    j = 0;
    nF = norm(f);
    tic;
    if solve_options.use_MATLAB
        L = ichol(A);
        R = L*L' - A;
        while norm(rk)>crit
            uk = L'\(L\(R*uk + f));
            rk = f - A*uk;
            j = j+1;
            ICBIM_conv(1,j) = norm(rk)/nF;
        end
    else
        L = IncompleteCholesky(A,n-1,dim);
        Dinv = inv(spdiags(spdiags(L,0),0,(n-1)^dim,(n-1)^dim));
        L1 = L*Dinv;
        R = L1*L' - A;
        while norm(rk)>crit
            uk = UpperSolver(L',LowerSolver(L1,R*uk + f,(n-1)^(dim-1)),(n-1)^(dim-1));
            rk = f - A*uk;
            j = j+1;
            ICBIM_conv(1,j) = norm(rk)/nF;
        end
    end
    times_IC = toc;


end

