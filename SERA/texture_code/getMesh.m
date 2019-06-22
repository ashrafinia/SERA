function [ROIboxF,ROIboxV] = getMesh(ROIbox,pixelWx,pixelWy,sliceS)
% This function gets an ROI of zeros and ones and generates triangular mesh
% surface based on Maching Cubes algorithm.

[dx,dy,dz]=size(ROIbox);
ROIp = padarray(ROIbox,[1 1 1]);

% added ceil to prevent dimension mistmach error for voxels with sub mm sizes
[x,y,z]=ndgrid(0:pixelWx:ceil(pixelWx*(dx+2)-1) , 0:pixelWy:ceil(pixelWy*(dy+2)-1) , 0:sliceS:ceil(sliceS*(dz+2)-1)); 
[ROIboxF,ROIboxV] = MarchingCubes(single(x),single(y),single(z),ROIp,0.5);