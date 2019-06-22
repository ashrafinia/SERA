function [GLRLM2D_Cmb, GLRLM2D_Avg] = getGLRLM3Dtex(ROIonly,levels)
% -------------------------------------------------------------------------
% function [GLRLM2D_Cmb, GLRLM2D_Avg] = getGLRLM3Dtex(ROIonly,levels)
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
% Note: ROIonly is the outputs of prepareVolume.m
% -------------------------------------------------------------------------
% OUTPUTS:
% - GLRLM3D_Cmb: 3D GLRLM features: First merging GLRLMs for all directions
%               then calculate features for the combined GLRLM matrix.
% - GLRLM3D_Avg: 3D GLRLM features calculate GLRLM features for each
%               direction, then average over all directions.
% -------------------------------------------------------------------------
% AUTHOR(S):    
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
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

%% Calculate GLRLM: 
[GLRLM, GLRLM_D] = getGLRLM(ROIonly , levels);


%% Calculating features 
GLRLM2D_Cmb = CalcGLRLM(GLRLM, ROIonly);        
FeatsTmp = CalcGLRLM(GLRLM_D, ROIonly);
GLRLM2D_Avg = mean(FeatsTmp , 2,'omitnan');



