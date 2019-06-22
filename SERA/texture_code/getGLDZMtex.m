function [GLDZM2D_F, GLDZM3D_F, GLDZM25D_F] = getGLDZMtex(ROI2D,ROI3D,ROIonly2D,ROIonly3D,levels2D,levels3D)
% -------------------------------------------------------------------------
% [GLDZM2D_F, GLDZM3D_F] =
% getGLDZMtex(ROI2D,ROI3D, ROIonly2D, ROIonly3D, levels2D, levels3D) 
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates the GLDZM matrix for 2D and 3D.
% In 2D, every slice is calculated separately, then features are calculated.
% 
% The grey level distance zone matrix (GLDZM) counts the number of groups
% of connected voxels with a specific discretised grey level value and
% distance to ROI edge (Thibault et al., 2014). The matrix captures the
% relation between location and grey level. Two maps are required to
% calculate the GLDZM. The first is a grey level grouping map, identical
% with the one created for the grey level size zone matrix (GLSZM). The
% second is a distance map.
% Note that for GLDZM we need both the intensity and morphological ROIs.
% -------------------------------------------------------------------------
% INPUTS:
% - ROI2D: Smallest box containing the 2D resampled ROI, with the imaging
%          data ready for texture analysis computations. Voxels outside the
%          ROI are set to NaNs.   
% - ROI3D: Smallest box containing the 3D resampled ROI, with the imaging
%          data ready for texture analysis computations. Voxels outside the
%          ROI are set to NaNs.   
% - ROIonly2D: Smallest box containing the 2D resampled morphological ROI,
%              values are either 0 or 1. 
% - ROIonly3D: Smallest box containing the 3D resampled morphological ROI,
%              values are either 0 or 1. 
% - levels2D: number of bins (for fixed number of bins method) or bin size
%           (for bin size method) for the 2D resampled ROI.
% - levels3D: number of bins (for fixed number of bins method) or bin size
%           (for bin size method) for the 3D resampled ROI.
% Note: ROIonly is the outputs of prepareVolume.m
% -------------------------------------------------------------------------
% OUTPUTS:
% - GLDZM2D: An array of 16 GLDZM features for the 2D resampled ROI.
% - GLDZM3D: An array of 16 GLDZM features for the 3D resampled ROI.
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

%% Calculate 2D GLDZM : 
[nX, nY, nZ] = size(ROI2D);
FeatTmp = [];
GLDZM_all = zeros(numel(levels2D) , ceil(max(nX,nY)/2) , nZ);
for s = 1:nZ
    [GLDZM]   = getGLDZM(ROI2D(:,:,s),ROIonly2D(:,:,s),levels2D);
    GLDZM_all(:,1:size(GLDZM,2) , s) = GLDZM; 
    
    % calculate 2D features
    GLDZMstr = CalcGLDZM(GLDZM,ROI2D(:,:,s));
    FeatTmp = cat(2,FeatTmp , GLDZMstr');
end

GLDZM2D_F = mean(FeatTmp , 2,'omitnan');


% Calculate 2.5D features
GLDZM25 = squeeze(sum(GLDZM_all, 3));
GLDZM25D_F = CalcGLDZM(GLDZM25,ROI2D)';


%% Calculate 3D GLDZM
GLDZM   = getGLDZM(ROI3D,ROIonly3D,levels3D);
[GLDZM3D_F] = CalcGLDZM(GLDZM,ROI3D)';


end