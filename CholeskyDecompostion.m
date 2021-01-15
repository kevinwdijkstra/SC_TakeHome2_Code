function [A] = CholeskyDecompostion(A)
% CHOLESYDECOMPOSITION Compute the Cholesky factor for matrix A according
% to Algorithm 1 from the lecture notes.
%   A = CHOLESKYDECOMPOSITION(A)  

    [n,~] = size(A); % Compute matrix dimensions n
    [p,~] = bandwidth(A); % Compute bandwidth of A
    for j = 1:n
        %% vectored
        Id1 = j+1:min(j+p,n);       % second for loop over i
        Id2 = max((j-p),1):(j-1);   % both sums over j
        
        % extract data in vector and matrix data so that all points over i
        % can be calculated at once (vectorized code for speed). This turns
        % the sum to a matrix vector product.
        Vec = A(j,Id2);
        Mat = A(Id1,Id2);
        
        % Note the norm is just a matrix vector product, however the norm
        % makes more clever use of the sparse vector for speed up.
        Cj     = sqrt(A(j,j)-norm(Vec)^2);
        A(j,j) = Cj; %store C(k,k) as accessing data takes some time.
        
        % calculate the entire second for loop over i in one go.
        A(Id1,j) = 1/Cj*(A(Id1,j) - Mat*Vec');
    end
    A = tril(A); % Extract lower triangular part of matrix

end

