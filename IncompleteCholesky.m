function [A] = IncompleteCholesky(A,p,dim)
% INCOMPLETECHOLESKY Compute Incomplete Cholesky factor with bands at
%   the diagonal and the off-diagonals +-1, +-p and optionally at +-q.

narginchk(2,3)
q = p^2; % Compute optional 3rd band
n = size(A,1); % Compute matrix dimension
A(n+1:(n+1):end) = 0; % Remove first band in the upper triangular part
A((n*p+1):(n+1):end) = 0; % Remove second band in the upper triangular part

for k=2:p % loop over diagonal elements until k = p
    A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1);
end

switch dim
    case 2 %2D case
        for k=(p+1):n % loop over diagonal elements from p to the end
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) - ...
                A(k,k-p)^2/A(k-p,k-p);
        end


    case 3 %3D case
        % Remove third band in the upper triangular part
        A((n*q+1):(n+1):end) = 0;
        
        for k=(p+1):q % loop over diagonal elements until k = p
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) ...
                - A(k,k-p)^2/A(k-p,k-p);
        end

        for k=(q+1):n % loop over diagonal elements from k = q to end
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) - ...
                A(k,k-p)^2/A(k-p,k-p) - A(k,k-q)^2/A(k-q,k-q);
        end
end

end