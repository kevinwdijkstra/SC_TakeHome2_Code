function [A] = IncompleteCholesky(A,p,q)
% INCOMPLETECHOLESKY Compute Incomplete Cholesky factor with bands at
%   the diagonal and the off-diagonals +-1, +-p and optionally at +-q.

narginchk(2,3)

n = size(A,1);
A(n+1:(n+1):end) = 0;
A((n*p+1):(n+1):end) = 0;

switch nargin
    case 2 %2D case
        for k=2:p
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1);
        end

        for k=(p+1):n
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) - ...
                A(k,k-p)^2/A(k-p,k-p);
        end


    case 3 %3D case
        A((n*q+1):(n+1):end) = 0;

        for k=(p+1):q
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) ...
                - A(k,k-p)^2/A(k-p,k-p);
        end

        for k=(q+1):n
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) - ...
                A(k,k-p)^2/A(k-p,k-p) - A(k,k-q)^2/A(k-q,k-q);
        end
end

end