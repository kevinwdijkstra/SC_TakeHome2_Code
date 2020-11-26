function [x] = LowerSolver(C,y)
%% solve Cx = y when C is a lower triangle matrix,

[~,m] = size(C);
x = zeros(1,m);


% [I,J,Val] = find(C);
for i = 1:m
%     Id  = (I==i);
%     j   = J(Id);
%     S   = Val(Id)';
    
    SubC = C(i,:);
    [~,j,S] = find(SubC);
    x(i) = (y(i) - sum(S(1:end-1).*x(j(1:end-1))))/S(end);
end



end

