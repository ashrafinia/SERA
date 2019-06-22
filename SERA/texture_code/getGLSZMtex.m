function [Feats_2D, Feats_3D, Feats_25D ] = getGLSZMtex(ROI2D,ROI3D,levels2D,levels3D)
% -------------------------------------------------------------------------
% [GLSZM2D, GLSZM3D] = getGLSZMtex(ROI2D,ROI3D,levels2D,levels3D)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates the GLSZM matrix for 2D and 3D.
% In 2D, every slice is calculated separately, then features are calculated.
% 
% The grey level size zone matrix (GLSZM) counts the number of groups of
% connected voxels witha specific discretised grey level value and size
% (Thibault et al., 2014). Voxels are connected ifthe neighbouring voxel
% has the same discretised grey level value.   
% -------------------------------------------------------------------------
% INPUTS:
% - ROI2D: Smallest box containing the 2D resampled ROI, with the imaging
%          data ready for texture analysis computations. Voxels outside the
%          ROI are set to NaNs.   
% - ROI3D: Smallest box containing the 3D resampled ROI, with the imaging
%          data ready for texture analysis computations. Voxels outside the
%          ROI are set to NaNs.   
% - levels2D: number of bins (for fixed number of bins method) or bin size
%           (for bin size method) for the 2D resampled ROI.
% - levels3D: number of bins (for fixed number of bins method) or bin size
%           (for bin size method) for the 3D resampled ROI.
% Note: ROIonly is the outputs of prepareVolume.m
% -------------------------------------------------------------------------
% OUTPUTS:
% - GLSZM2D: An array of 16 GLSZM features for the 2D resampled ROI.
% - GLSZM3D: An array of 16 GLSZM features for the 3D resampled ROI.
% -------------------------------------------------------------------------
% AUTHOR(S):    
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2017
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of Radiomics Package by Saeed Ashrafinia, Rahmimlab.com
% --> Copyright (C) 2013-2017  Saeed Ashrafinia, Johns Hopkins University
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


%% Calculate 2D GLSZM : 
[nX, nY, nZ] = size(ROI2D);
FeatTmp = [];
GLSZM2D_all = zeros(numel(levels2D) , ceil(max(nX,nY)/2) , nZ); 
for s = 1:nZ
    [GLSZM]   = getGLSZM(ROI2D(:,:,s),levels2D);
    GLSZM2D_all(:,1:size(GLSZM,2),s) = GLSZM;
    
    % calculate 2D feature on-the-fly
    GLSZMstr = struct2array(getGLSZMtextures(GLSZM));
    FeatTmp = cat(2,FeatTmp , GLSZMstr');
    
end

Feats_2D = mean(FeatTmp , 2,'omitnan');

% Calculate 2.5D features
GLSZM25 = sum(GLSZM2D_all, 3);
if sum(GLSZM25(:))~=0
    Feats_25D = struct2array(getGLSZMtextures(GLSZM25))';
else
    Feats_25D = Feats_2D;
end



%% Calculate 3D GLSZM
[GLSZM]   = getGLSZM(ROI3D,levels3D);
[GLSZM3Dstr] = getGLSZMtextures(GLSZM);
Feats_3D = struct2array(GLSZM3Dstr)';


end