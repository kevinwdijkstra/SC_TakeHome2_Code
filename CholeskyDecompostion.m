function [C] = CholeskyDecompostion(A)
   % Algotihm 4 from book


    [n,~] = size(A);
%     [i,j,a] = find(A);
    [p,~] = bandwidth(A);
    
    C = 0*A;
    for j = 1:n
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
        %% vectored
        Id1 = j+1:min(j+p,n);       % second for loop over i
        Id2 = max((j-p),1):(j-1);   % both sums over j
        
        % extract data in vector and matrix data so that all points over i
        % can be calculated at once (vectorized code for speed). This turns
        % the sum to a matrix vector product.
        Vec = C(j,Id2);
        Mat = C(Id1,Id2);
        
        % Note the norm is just a matrix vector product, however the norm
        % makes more clever use of the sparse vector for speed up.
        Cj     = sqrt(A(j,j)-norm(Vec)^2);
        C(j,j) = Cj; %store C(k,k) as accessing data takes some time.
        
        % calculate the entire second for loop over i in one go.
        C(Id1,j) = 1/Cj*(A(Id1,j) - Mat*Vec');
    end

end

