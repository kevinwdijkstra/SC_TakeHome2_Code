function [y] = LowerSolver(C,y,p)
%% solve Cx = y when C is a lower triangle matrix,
% inputs:
% - C is the lower triangle matrix
% - y is the right hand side vector
% - p is the band width of the matrix C
% we reuse the y vector to store the solution to save on memory.

% check amount of inputs
narginchk(2,3)

% find size of matrix C
[~,m] = size(C);

% we use a different algorithm when the bandwidth is known
switch nargin
    case 2 % we DONT know the bandwidth
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
        
    case 3 % we DO know the bandwidth
        % loop over y from first value to last.
        for i = 1:m
            % calculate over how many values we have to loop
            i_min = max(1,i-p);
            % get the ith row from C
            SubC = C(i,i_min:i);
            % calculate new C. This is x(i) = (y(i) - sum_j=1^i-1 C(i,j)x(j))/C(i,i)
            % By indexing only the necesary values we get a speed up
            y(i) = (y(i) - sum(SubC(1:(end-1))*y(i_min:(i-1))))/SubC(end);
        end      
end

end

