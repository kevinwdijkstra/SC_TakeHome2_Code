function [x] = LowerSolver(C,y)
%% solve Cx = y when C is a lower triangle matrix,
% find size of matrix C and pre-alocate memory
[~,m] = size(C);
x = zeros(1,m);

% loop over x from first value to last.
for i = 1:m
    % get the ith row from C
    SubC = C(i,:);
    
    % get the y index and value of the ith row of C
    [~,j,S] = find(SubC);
    
    % calculate new C. This is x(i) = (y(i) - sum_j=1^i-1 C(i,j)x(j))/C(i,i)
    % By indexing only the necesary values we get a speed up
    x(i) = (y(i) - sum(S(1:end-1).*x(j(1:end-1))))/S(end);
end



end

