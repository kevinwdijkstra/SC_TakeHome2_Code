function [D2Mat] = CreateMatrix2D(n)

    D1Mat = Create1D(n);

    D2Mat = (kron(D1Mat,speye(n-1))+kron(speye(n-1),D1Mat));
    
    D2Mat = (n)^2*D2Mat;

    
end


function [D1Mat] = Create1D(n)
    e = ones(n-1,1);
    D1Mat = spdiags([-e 2*e -e],-1:1,n-1,n-1);
end

