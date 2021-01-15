function [A] = CholeskyDecompostion(A)
%% calculate the Cholesky Decomposition
% inputs:
% - A the matrix

% This code is an implementation of algotihm 4 of the  book

% get matrix size
[n,~] = size(A);
% get bandwidth
[p,~] = bandwidth(A);

for j = 1:n % loop over  each row
    %% the code we use is the following sequential code vectorized for MATLAB speed up
%         Ck     = sqrt(A(k,k)-sum(C(k,max((k-N),1):k-1).^2));
%         C(k,k) = Ck;
%         A(k,k) = Ck;
%         for i = k+1:min(k+M,n)
%             
%             Csub   = sum(C(i,max((k-N),1):(k-1)).*C(k,max((k-N),1):(k-1)));
%             C(i,k) = 1/Ck*(A(i,k)-Csub);
%             
%             
%             
%             A(i,k) = C(i,k);
%         end

    %% vectored implementation for speed up:
    Id1 = j+1:min(j+p,n);       % second for loop over i in line 18
    Id2 = max((j-p),1):(j-1);   % both sums over j in line 20

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
A = tril(A);

end

