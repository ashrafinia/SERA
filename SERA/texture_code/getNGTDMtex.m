function [NGTDM2D, NGTDM3D, NGTDM25D] = getNGTDMtex(ROI2D,ROI3D,levels2D,levels3D)
% -------------------------------------------------------------------------
% [NGTDM2D, NGTDM3D] = getNGTDMtex(ROI2D,ROI3D,levels2D,levels3D)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates the NGTDM matrix for 2D and 3D.
% In 2D, every slice is calculated separately, then features are calculated.
% 
% The grey level size zone matrix (NGTDM) contains the sum of grey level
% differences of pixels/voxels with discretised grey level i and the
% average discretised grey level of neighbouring pixels/voxels within a
% distance d. 
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
% - NGTDM2D: An array of 16 NGTDM features for the 2D resampled ROI.
% - NGTDM3D: An array of 16 NGTDM features for the 3D resampled ROI.
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


%% Calculate 2D NGTDM : 
[~, ~, nZ] = size(ROI2D);
FeatTmp = [];
NGTDM2D_all = zeros(numel(levels2D) , nZ); count_valid_all = NGTDM2D_all;

for s = 1:nZ
    [NGTDM,countValid]   = getNGTDM(ROI2D(:,:,s),levels2D);
    NGTDM2D_all(:,s) = NGTDM;
    count_valid_all(:,s) = countValid; 
    
    % calculate 2D feature on-the-fly
    NGTDMstr = struct2array(getNGTDMtextures(NGTDM,countValid));
    FeatTmp = cat(2,FeatTmp , NGTDMstr');
end

NGTDM2D = mean(FeatTmp , 2,'omitnan');

% Calculate 2.5D features
NGTDM25 = sum(NGTDM2D_all, 2);
NGTDM25D = struct2array(getNGTDMtextures(NGTDM25,sum(count_valid_all,2)))';


%% Calculate 3D NGTDM
[NGTDM,countValid]   = getNGTDM(ROI3D,levels3D);
[NGTDM3Dstr] = getNGTDMtextures(NGTDM,countValid);
NGTDM3D = struct2array(NGTDM3Dstr)';


end