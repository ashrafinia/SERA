function [area_out] = getSurfaceMesh(F,V)

% This function calculates the surface area of a mesh created with Maching
% Cubes algorithm with F and V.

a = V(F(:, 2), :) - V(F(:, 1), :);
b = V(F(:, 3), :) - V(F(:, 1), :);
c = cross(a, b, 2);
area_out = 1/2 * sum(sqrt(sum(c.^2, 2)));