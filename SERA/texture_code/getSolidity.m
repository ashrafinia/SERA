function [solidity , CHarea] = getSolidity(ROIonly,pixelW,sliceS)
% -------------------------------------------------------------------------
% function [solidity] = getSolidity(ROIonly,pixelW,sliceS)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the solidity metric of the region of interest 
% (ROI) of an input volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% - pixelW: Pixel width, or in-plane resolution, in mm.
% - sliceS: Slice spacing, in mm.
% -------------------------------------------------------------------------
% OUTPUTS:
% - solidity: Ratio of the number of voxels in the ROI to the number of 
%             voxels in the 3D convex hull of the ROI (smallest polyhedron 
%             containing the ROI).
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres
%
%    This package is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this package.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

mask = ~isnan(ROIonly); % Find mask covering the ROI

% ISOTROPIC RESAMPLING
sFactor = sliceS/pixelW; % scaling factor
% mask = imresize3D(mask,[],[round(double(size(mask,1))),round(double(size(mask,2))),round(double(size(mask,3))*sFactor)],'nearest','fill');

% Replacing the above line with this one to take care of errors happenning
% e.g. for [1 0; 0 1]
try
    mask = imresize3D(mask,[],[round(double(size(mask,1))),round(double(size(mask,2))),round(double(size(mask,3))*sFactor)],'nearest','fill');
catch
    warning(['Solidity: adding a third dimension to "mask" with size ',mat2str(size(mask)), ' and values of  ',mat2str(mask(:))]);
    mask = imresize3D(cat(3,mask,NaN(size(mask))),[],[round(double(size(mask,1))),round(double(size(mask,2))),round(double(size(mask,3))*sFactor)],'nearest','fill');
end

% SOLIDITY COMPUTATION
perimeter = bwperim(mask,18);
nPoints = length(find(perimeter));
X = zeros(nPoints,1); Y = zeros(nPoints,1); Z = zeros(nPoints,1);
count = 1;
for i = 1:size(perimeter,3)
    [row,col] = find(perimeter(:,:,i));
    p = length(row);
    if p > 0
        X(count:count+p-1,1) = col(1:end);
        Y(count:count+p-1,1) = row(1:end);
        Z(count:count+p-1,1) = i;
        count = count + p;
    end
end
if peak2peak(Z)>1 && peak2peak(Y)>1 && peak2peak(X)>1
    [K,volumeH] = convhull(X,Y,Z);
else
    try
        [K,volumeH] = convhull(X,Y);
    catch
        volumeH = 0; 
    end
    
end

    
volumeROI = sum(mask(:));
solidity = volumeROI/volumeH;
if isinf(solidity), solidity = 0; end 

% Take care of K when it has just one dimension:
try
    if numel(K)==max(size(K))
        K(:,2:3) = ones(size(K,1),2);
    end
catch
    CHarea = 0;
    return;
end

CHarea = sum(sqrt(sum(( ...
 [Y(K(:,1)).*Z(K(:,2)) - Z(K(:,1)).*Y(K(:,2)) ...
  Z(K(:,1)).*X(K(:,2)) - X(K(:,1)).*Z(K(:,2))  ...
  X(K(:,1)).*Y(K(:,2)) - Y(K(:,1)).*X(K(:,2))] + ...
 [Y(K(:,2)).*Z(K(:,3)) - Z(K(:,2)).*Y(K(:,3)) ...
  Z(K(:,2)).*X(K(:,3)) - X(K(:,2)).*Z(K(:,3))  ...
  X(K(:,2)).*Y(K(:,3)) - Y(K(:,2)).*X(K(:,3))] + ...
 [Y(K(:,3)).*Z(K(:,1)) - Z(K(:,3)).*Y(K(:,1)) ...
  Z(K(:,3)).*X(K(:,1)) - X(K(:,3)).*Z(K(:,1))  ...
  X(K(:,3)).*Y(K(:,1)) - Y(K(:,3)).*X(K(:,1))]).^2,2))) / 2;

end