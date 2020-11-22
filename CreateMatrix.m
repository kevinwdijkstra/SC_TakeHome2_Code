function [D2Mat, D3Mat] = CreateMatrix(N)

    D1Mat = Create1D(N);

    D2Mat = kron(D1Mat,speye(N-1))+kron(speye(N-1),D1Mat);
    D3Mat = kron(D2Mat,speye(N-1))+kron(speye((N-1)^2),D1Mat);
    
    D2Mat = (N)^2*D2Mat;
    D3Mat = (N)^3*D3Mat;
    
end




function [D1Mat] = Create1D(N)
    e = ones(N-1,1);
    D1Mat = spdiags([-e 2*e -e],-1:1,N-1,N-1);
end

