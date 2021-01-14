function [D2f_dir] = CreateBC2D(u,f,Mesh,n)
%% calculate the right hand side vector for the poisson problem with boundary elimination
% inputs:
% - u the exact solution function
% - f the source function
% - Mesh the mesh
% - n the number of elements in each dimension.


% calculate the source term.
D2f = f(Mesh);

% calculate logicals for which elements are close to the boundary
Ixlow = (Mesh(1,:)-1/(n)==0);
Ixhig = (Mesh(1,:)+1/(n)==1);
Iylow = (Mesh(2,:)-1/(n)==0);
Iyhig = (Mesh(2,:)+1/(n)==1);

% encorporate the boundary values according to the logicals
D2f_dir     = D2f + n^2*(   Ixlow.*u(Mesh - [1/(n);0])+...
                            Ixhig.*u(Mesh + [1/(n);0])+...
                            Iylow.*u(Mesh - [0;1/(n)])+...
                            Iyhig.*u(Mesh + [0;1/(n)]));
end

