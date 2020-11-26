function [D2Mesh, x, y, z] = CreateMesh2D(n)

    x = linspace(0,1,n+1);
    y = linspace(0,1,n+1);
    z = linspace(0,1,n+1);
    
    x = x(2:end-1);
    y = y(2:end-1);
    z = z(2:end-1);
    
    D2Mesh = combvec(x,y);
end

