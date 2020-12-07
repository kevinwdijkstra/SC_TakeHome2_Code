function [A] = IncompleteCholesky(A,p)
    n = size(A,1);
    A(n+1:(n+1):end) = 0;
    A((n*p(1)+1):(n+1):end) = 0;
    if numel(p)== 2
        A((n*p(2)+1):(n+1):end) = 0;
    end
    
    
    for k=2:p(1)
        A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1);
    end
    
    
    if numel(p) == 1
        for k=(p(1)+1):n
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) - A(k,k-p(1))^2/A(k-p(1),k-p(1));
        end
        
    else
        for k=(p(1)+1):p(2)
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) - A(k,k-p(1))^2/A(k-p(1),k-p(1));
        end
        
        for k=(p(2)+1):n
            A(k,k) = A(k,k) - A(k,k-1)^2/A(k-1,k-1) - A(k,k-p(1))^2/A(k-p(1),k-p(1)) - A(k,k-p(2))^2/A(k-p(2),k-p(2));
        end
    end
    
    
    
end