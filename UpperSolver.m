function [x] = UpperSolver(C,y)
%% solve Cx = y when C is a upper triangle matrix,

[~,m] = size(C);
x = zeros(1,m);

for i = m:-1:1
    SubC = C(i,:);
    [~,j,S] = find(SubC);
    x(i) = (y(i) - sum(S(2:end).*x(j(2:end))))/S(1);
end



end

