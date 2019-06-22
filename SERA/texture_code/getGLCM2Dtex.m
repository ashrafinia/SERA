function [Feats_KSKD, Feats_MSKD, Feats_KSMD, Feats_MSMD] = getGLCM2Dtex(ROIonly,levels)
% -------------------------------------------------------------------------
% function [SM_f, SS_f] = getGLCM2Dtex(ROIonly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates 2D GLCM and calculates texture features. 
% The beginning of this code is from Martin Valleries GLCM code 
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
% - SM_f: 2D GLCM features merging GLCMs of slices, then calculate features
% - SS_f: 2D GLCM features calculate GLCM features for each slice, then
%         average over slices
% -------------------------------------------------------------------------
% AUTHOR(S):    
% - Saeed Ashrafinia
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013.
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
%% Initialization
nLevel = length(levels);
if nLevel > 100, adjust = 10000; else, adjust = 1000; end
levelTemp = max(levels)+1;
ROIonly(isnan(ROIonly)) = levelTemp;
levels = [levels,levelTemp];

dim = size(ROIonly);
if ndims(ROIonly) == 2 %#ok<ISMAT>
	dim(3) = 1;
end
q2 = reshape(ROIonly,1,prod(dim));


% QUANTIZATION EFFECTS CORRECTION (M. Vallieres)
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
qs = round(levels*adjust)/adjust;
q2 = round(q2*adjust)/adjust;

q3 = q2*0;
for k = 1:length(qs)
	q3(q2==qs(k)) = k;
end
ROInanReplaced = reshape(q3,dim);

%% Create GLCMs
DIST = [[1 0];[1 1];[0 1];[-1 1]];
% DIST = [[0 1];[1 0];[1 1];[0 -1];[0 -1];[-1 0];[-1 1];[1 -1]];

[~, ~, nZ] = size(ROInanReplaced);
nlevels = length(levels)-1 ;% minus one to only save the GLCM for levels included in ROI. It should only be +1 when calling graycomatrix function below.
GLCM_KeepSlice_KeepDirs = zeros(nlevels,nlevels,size(DIST,1),nZ);
if ~isempty(find(ROInanReplaced==0)), error('Need to change the minimum GrayLimits of GLCM to zeros.'); end %#ok<EFIND>
for s = 1:nZ
    tmpGLCM = graycomatrix(ROInanReplaced(:,:,s),'Offset',DIST,'NumLevels',nlevels+1,'GrayLimits',[1 nlevels+1],'Symmetric', true);
    tmpGLCM = tmpGLCM(1:(end-1) , 1:(end-1),:);
    GLCM_KeepSlice_KeepDirs(:,:,:,s) =  tmpGLCM; 
    try
        if s ==1, warning('off','last'); end
    catch
    end
end

%% Create GLCMs based on keeping/merging slices and directions
% Merge Slices, Merge Dirs
GLCM_AllMerged   = sum(sum(GLCM_KeepSlice_KeepDirs,3),4);
GLCMnorm_AllMerged = GLCM_AllMerged / sum(GLCM_AllMerged(:)); % Normalize GLCM

% Keep Slices, Merge Dirs
GLCM_KeepSlice_MergeDirs = zeros(size(GLCM_AllMerged,1), size(GLCM_AllMerged,2),size(GLCM_KeepSlice_KeepDirs,4));
GLCM_KeepSlice_MergeDirs(:,:,:) = sum(GLCM_KeepSlice_KeepDirs , 3);

tmp = sum(sum(GLCM_KeepSlice_MergeDirs , 1),2);
GLCMnorm_KeepSlice_MergeDirs = GLCM_KeepSlice_MergeDirs ./ repmat(tmp ,size(GLCM_KeepSlice_MergeDirs,1),size(GLCM_KeepSlice_MergeDirs,2),1);
GLCMnorm_KeepSlice_MergeDirs(isnan(GLCMnorm_KeepSlice_MergeDirs)) = 0;

% Keep Slices, Keep Dirs
tmp = sum(sum(GLCM_KeepSlice_KeepDirs , 1),2);
GLCMnorm_KeepSlice_KeepDirs = GLCM_KeepSlice_KeepDirs ./ repmat(tmp ,size(GLCM_KeepSlice_MergeDirs,1),size(GLCM_KeepSlice_MergeDirs,2),1,1);
GLCMnorm_KeepSlice_KeepDirs(isnan(GLCMnorm_KeepSlice_KeepDirs)) = 0;

% Merge Slices, Keep Dirs
GLCM_MergeSlice_KeepDirs = zeros(size(GLCM_AllMerged,1), size(GLCM_AllMerged,2),size(GLCM_KeepSlice_KeepDirs,3));
GLCM_MergeSlice_KeepDirs(:,:,:) = sum(GLCM_KeepSlice_KeepDirs , 4);

tmp = sum(sum(GLCM_MergeSlice_KeepDirs , 1),2);
GLCMnorm_MergeSlice_KeepDirs = GLCM_MergeSlice_KeepDirs ./ repmat(tmp ,size(GLCM_KeepSlice_KeepDirs,1),size(GLCM_KeepSlice_KeepDirs,2),1);
GLCMnorm_MergeSlice_KeepDirs(isnan(GLCMnorm_MergeSlice_KeepDirs)) = 0;


%% Calculating features 

% Calculate features for AllMerged
Feats_MSMD = CalcGLCM(GLCMnorm_AllMerged);
nFeats = size(Feats_MSMD,1);

% Calculate features for Keep Slices, Merge Dirs
tmp = CalcGLCM(GLCMnorm_KeepSlice_MergeDirs);
if ndims(tmp)>2, warning('GLCM features might having an extra dimension'); end %#ok<ISMAT>
Feats_KSMD = mean(tmp , 2,'omitnan');

% Calculate features for Merge Slices, Keep Dirs
tmp = CalcGLCM(GLCMnorm_MergeSlice_KeepDirs);
if ndims(tmp)>2, warning('GLCM features might having an extra dimension'); end %#ok<ISMAT>
Feats_MSKD = mean(tmp , 2,'omitnan');

% Calculate features for Keep Slices, Keep Dirs
tmp = zeros(nFeats , size(DIST,1), nZ);
for d = 1:size(DIST,1)
    for s = 1:nZ
        tmp(:,d,s) = CalcGLCM(squeeze(GLCMnorm_KeepSlice_KeepDirs(:,:,d,s)));
    end
end
Feats_KSKD = mean(squeeze(mean(tmp , 3,'omitnan')),2,'omitnan');



end






