function [C] = CholeskyDecompostion(A)
   % Algotihm 4 from book


    [n,~] = size(A);
%     [i,j,a] = find(A);
    [N,M] = bandwidth(A);
    
    C = 0*A;
    for k = 1:n
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
        Id1 = k+1:min(k+M,n);
        Id2 = max((k-N),1):(k-1);
        
        Vec = C(k,Id2);
        Mat = C(Id1,Id2);
        
        
        Ck     = sqrt(A(k,k)-norm(Vec)^2);
        C(k,k) = Ck;
        
        C(Id1,k) = 1/Ck*(A(Id1,k) - Mat*Vec');
%         A(Id1,k) = C(Id1,k);
        
        
    end

end

