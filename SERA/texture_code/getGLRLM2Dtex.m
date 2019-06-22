function [Feats_KSKD, Feats_MSKD, Feats_KSMD, Feats_MSMD] = getGLRLM2Dtex(ROIonly,levels)
% -------------------------------------------------------------------------
% function [GLRLM2D_Cmb, GLRLM2D_Avg] = getGLRLM2Dtex(ROIonly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates the GLRLM matrix for each slice of an ROI.
% Every slice is calculated separately, then features are calculated.
% 
% Like the grey level co-occurrence matrix, GLRLM also assesses the
% distribution of discretised grey levels in an image or in a stack of
% images. However, instead of assessing the combination of levels between
% neighbouring pixels or voxels, GLRLM assesses grey level run % lengths.
% Run length counts the frequency of consecutive voxels with discretised
% grey level i along  direction \Delta.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: Smallest box containing the ROI, with the imaging data ready 
%            for texture analysis computations. Voxels outside the ROI are 
%            set to NaNs. 
% - levels: number of bins (for fixed number of bins method) or bin size
%           (for bin size method) 
% Note: RIOonly is the outputs of prepareVolume.m
% -------------------------------------------------------------------------
% OUTPUTS:
% - GLRLM2D_Cmb: 2D GLRLM features: merging GLRLM for different slice of
%                 the volume, then calculate features for the combined
%                 GLRLM matrix. 
% - GLRLM2D_Avg: 2D GLRLM features: calculate GLRLM features for each
%                slice first, then average over all slices.
% -------------------------------------------------------------------------
% AUTHOR(S):    
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
% - Revision: July 2017
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

%% Calculate GLCM: 
% The following function returns a 3D matrix of GLCMs for all 13 directions
[nx, ny, nZ] = size(ROIonly);
FeatTmp = [];
FeatTmp_D = [];
GLRLM_all = zeros(numel(levels) , max(nx,ny) , nZ);
GLRLM_D_all = zeros(numel(levels) , max(nx,ny) , 4 , nZ);
for s = 1:nZ
    [GLRLM, GLRLM_D] = getGLRLM(ROIonly(:,:,s) , levels);
    FeatTmp = cat(2,FeatTmp , CalcGLRLM(GLRLM, ROIonly(:,:,s)));
    FeatTmp_D = cat(3,FeatTmp_D , CalcGLRLM(GLRLM_D(:,:,1:4), ROIonly(:,:,s)));
    GLRLM_all(:,1:size(GLRLM,2),s) = GLRLM;
    GLRLM_D_all(:,1:size(GLRLM_D,2),1:4,s) = GLRLM_D(:,:,1:4);
    
end


% Keep Slices, Merge Dirs
Feats_MSKD = mean(FeatTmp , 2,'omitnan');
Feats_MSKD(13) = Feats_MSKD(13) /4;

% Keep slice, Keep Dirs
Feats_KSKD = mean(reshape(FeatTmp_D,16,[]),2,'omitnan');

% Merge Slices, Keep Dirs
GLRLM_MergeSlice_KeepDirs = zeros(size(GLRLM_D_all,1), size(GLRLM_D_all,2),size(GLRLM_D_all,3));
GLRLM_MergeSlice_KeepDirs(:,:,:) = sum(GLRLM_D_all , 4);

tmp = sum(sum(GLRLM_MergeSlice_KeepDirs , 1),2);
GLRLMnorm_MergeSlice_KeepDirs = GLRLM_MergeSlice_KeepDirs ./ repmat(tmp ,size(GLRLM_D_all,1),size(GLRLM_D_all,2),1);
GLRLMnorm_MergeSlice_KeepDirs(isnan(GLRLMnorm_MergeSlice_KeepDirs)) = 0;

Feats_KSMD = CalcGLRLM(GLRLMnorm_MergeSlice_KeepDirs, ROIonly);
Feats_KSMD = mean(Feats_KSMD , 2,'omitnan');


% Merge Slices, Merge Dirs
GLRLM_AllMerged   = sum(sum(GLRLM_D_all,3),4);
GLRLMnorm_AllMerged = GLRLM_AllMerged / sum(GLRLM_AllMerged(:)); % Normalize GLCM

Feats_MSMD = CalcGLRLM(GLRLMnorm_AllMerged, ROIonly);
Feats_MSMD = mean(Feats_MSMD , 2,'omitnan');



end