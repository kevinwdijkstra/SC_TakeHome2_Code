function [D3f_dir] = CreateBC3D(u,f,Mesh,n)
D3f = f(Mesh);

Ixlow = (Mesh(1,:)-1/(n)==0);
Ixhig = (Mesh(1,:)+1/(n)==1);
Iylow = (Mesh(2,:)-1/(n)==0);
Iyhig = (Mesh(2,:)+1/(n)==1);
Izlow = (Mesh(3,:)-1/(n)==0);
Izhig = (Mesh(3,:)+1/(n)==1);

D3f_dir     = D3f + n^2*(   Ixlow.*u(Mesh - [1/(n);0;0])+...
                            Ixhig.*u(Mesh + [1/(n);0;0])+...
                            Iylow.*u(Mesh - [0;1/(n);0])+...
                            Iyhig.*u(Mesh + [0;1/(n);0])+...
                            Izlow.*u(Mesh - [0;0;1/(n)])+...
                            Izhig.*u(Mesh + [0;0;1/(n)]));
end

