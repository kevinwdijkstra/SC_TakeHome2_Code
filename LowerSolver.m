function [y] = LowerSolver(C,y,p)
narginchk(2,3)
%% solve Cx = y when C is a lower triangle matrix,
% find size of matrix C and pre-alocate memory
[~,m] = size(C);
%x = zeros(1,m);

switch nargin
    case 2
        % loop over x from first value to last.
        for i = 1:m
            % get the ith row from C
            SubC = C(i,1:i);

            % get the y index and value of the ith row of C
            [~,j,S] = find(SubC);

            % calculate new C. This is x(i) = (y(i) - sum_j=1^i-1 C(i,j)x(j))/C(i,i)
            % By indexing only the necesary values we get a speed up
            y(i) = (y(i) - sum(S(1:end-1)*y(j(1:end-1))))/S(end);
        end
        
    case 3
        % loop over y from first value to last.
        for i = 1:m
            i_min = max(1,i-p);
            SubC = C(i,i_min:i);
            y(i) = (y(i) - sum(SubC(1:(end-1))*y(i_min:(i-1))))/SubC(end);
        end
        
        
end



end

