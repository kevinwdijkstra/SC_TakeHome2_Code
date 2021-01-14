function [D3Mesh, x, y, z] = CreateMesh3D(n)
%% calculate the mesh and x,y,z values of the unit square/cube with n 
% elements in each dimension
% inputs:
% - n the number of elements in each dimension.

    % calculate the x,y,z coordinates
    x = linspace(0,1,n+1);
    y = linspace(0,1,n+1);
    z = linspace(0,1,n+1);
    
    % remove the boundary values
    x = x(2:end-1);
    y = y(2:end-1);
    z = z(2:end-1);
    
    % create a 2xN vector of all possible position combinations.
    D2Mesh = combvec(x,y);
    % create a 3xN vector of all possible position combinations.
    D3Mesh = combvec(D2Mesh,z);


end

