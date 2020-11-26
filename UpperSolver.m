function [x] = UpperSolver(C,y)
%% solve Cx = y when C is a upper triangle matrix,
% find size of matrix C and pre-alocate memory
[~,m] = size(C);
x = zeros(1,m);

% loop over x from last value to first.
for i = m:-1:1
    % get the ith row from C
    SubC = C(i,:);
    
    % get the y index and value of the ith row of C
    [~,j,S] = find(SubC);
    
    % calculate new C. This is x(i) = (y(i) - sum_j=1^i-1 C(i,j)x(j))/C(i,i)
    % By indexing only the necesary values we get a speed up
    x(i) = (y(i) - sum(S(2:end).*x(j(2:end))))/S(1);
end



end

