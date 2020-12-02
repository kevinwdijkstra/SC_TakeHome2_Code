function [D2f_dir] = CreateBC2D(u,f,Mesh,n)
D2f = f(Mesh);

Ixlow = (Mesh(1,:)-1/(n)==0);
Ixhig = (Mesh(1,:)+1/(n)==1);
Iylow = (Mesh(2,:)-1/(n)==0);
Iyhig = (Mesh(2,:)+1/(n)==1);

D2f_dir     = D2f + n^2*(   Ixlow.*u(Mesh - [1/(n);0])+...
                            Ixhig.*u(Mesh + [1/(n);0])+...
                            Iylow.*u(Mesh - [0;1/(n)])+...
                            Iyhig.*u(Mesh + [0;1/(n)]));
end

