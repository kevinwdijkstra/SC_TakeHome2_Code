function [y] = UpperSolver(C,y,p)
%% solve Cx = y when C is a upper triangle matrix,
% inputs:
% - C is the upper triangle matrix
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
        
    case 3 % we DO know the bandwidth
        % loop over x from last value to first.
        for i = m:-1:1
            % calculate over how many values we have to loop
            i_max = min(m,i+p);
            % get the ith row from C
            SubC = C(i,i:i_max);
            % calculate new C. This is x(i) = (y(i) - sum_j=1^i-1 C(i,j)x(j))/C(i,i)
            % By indexing only the necesary values we get a speed up
            y(i) = (y(i) - sum(SubC(2:end)*y((i+1):i_max)))/SubC(1);
        end  
end

end

