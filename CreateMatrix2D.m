function [D2Mat] = CreateMatrix2D(n)
%% calculate the Matrix for the Poisson problem with boundary elimination.
% inputs:
% - n the number of elements in each dimension.

    % calculate the 1D poisson matrix
    D1Mat = Create1D(n);
    
    % calculate the 2D poisson matrix from the 1D with tensor products
    D2Mat = (kron(D1Mat,speye(n-1))+kron(speye(n-1),D1Mat));
    
    % scale the matrix correctly
    D2Mat = (n)^2*D2Mat;
end


function [D1Mat] = Create1D(n)
    % set up a ones vector
    e = ones(n-1,1);
    % create a 1D poisson sparse matrix
    D1Mat = spdiags([-e 2*e -e],-1:1,n-1,n-1);
end

