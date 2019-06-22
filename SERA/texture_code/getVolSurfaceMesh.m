function [vol_out, area_out] = getVolSurfaceMesh(F,V)

% This function calculates the volume and surface area of a mesh created
% with Maching Cubes algorithm with F and V.

%% Calculate volume from mesh
% initialize an array of volume
nFaces = size(F, 1);
vols = zeros(nFaces, 1);

% Shift all vertices to the mesh centroid
vertices = bsxfun(@minus, V, mean(V,1));

% compute volume of each tetraedron
for i = 1:nFaces
    % consider the tetrahedron formed by face and mesh centroid
    tetra = vertices(F(i, :), :);
    
    % volume of current tetrahedron
    vols(i) = det(tetra) / 6;
end

vol_out = abs(sum(vols));


%% Calculate the surface
a = V(F(:, 2), :) - V(F(:, 1), :);
b = V(F(:, 3), :) - V(F(:, 1), :);
c = cross(a, b, 2);
area_out = 1/2 * sum(sqrt(sum(c.^2, 2)));