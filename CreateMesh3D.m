function [D3Mesh, x, y, z] = CreateMesh3D(n)

    x = linspace(0,1,n+1);
    y = linspace(0,1,n+1);
    z = linspace(0,1,n+1);
    
    x = x(2:end-1);
    y = y(2:end-1);
    z = z(2:end-1);
    
    D2Mesh = combvec(x,y);
    D3Mesh = combvec(D2Mesh,z);


end

