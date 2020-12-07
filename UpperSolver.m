function [y] = UpperSolver(C,y,p)
narginchk(2,3)
%% solve Cx = y when C is a upper triangle matrix,
% find size of matrix C and pre-alocate memory
[~,m] = size(C);
%x = zeros(1,m);

switch nargin
    case 2
        % loop over x from last value to first.
        for i = m:-1:1
            % get the ith row from C
            SubC = C(i,i:end);

            % get the y index and value of the ith row of C
            [~,j,S] = find(SubC);

            % calculate new C. This is x(i) = (y(i) - sum_j=1^i-1 C(i,j)x(j))/C(i,i)
            % By indexing only the necesary values we get a speed up
            y(i) = (y(i) - sum(S(2:end)*y(j(2:end))))/S(1);
        end
        
    case 3
        for i = m:-1:1
            i_max = min(m,i+p);
            SubC = C(i,i:i_max);
            y(i) = (y(i) - sum(SubC(2:end)*y((i+1):i_max)))/SubC(1);
        end
        
        
end



end

