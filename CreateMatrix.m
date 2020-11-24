function [D2Mat, D3Mat] = CreateMatrix(n)

    D1Mat = Create1D(n-1);

    D2Mat = (kron(D1Mat,eye(n-1))+kron(eye(n-1),D1Mat));
    D3Mat = kron(D2Mat,eye(n-1))+kron(eye((n-1)^2),D1Mat);
    
    D2Mat = (n)^2*D2Mat;
    D3Mat = (n)^3*D3Mat;
    
end




function [D1Mat] = Create1D(n)
%     D1Mat = (2*eye(n) -diag(ones(n-1,1),1) - diag(ones(n-1,1),-1));
    D1Mat = spdiags([-1 2 -1].*ones(n,1),-1:1,n,n);
end

