function [ROIonlyUniform,levels] = fixedBinSizeQuantization(ROIonly,BinSize,minGL)
% -------------------------------------------------------------------------
% function [ROIonlyUniform,levels] = uniformQuantization(ROIonly,Ng)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes uniform quantization on the region of interest 
% (ROI) of an input volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% - Ng: Integer specifying the number of gray levels in the quantization.
% -------------------------------------------------------------------------
% OUTPUTS:
% - ROIonlyUniform: Quantized input volume.
% - levels: Vector containing the quantized gray-level values in the ROI
%           (or reconstruction levels of quantization).
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: Dec 2017
% -------------------------------------------------------------------------
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

ROIonly = double(ROIonly);
% maxVal = max(ROIonly(:));
% minVal = floor(min(ROIonly(:)));
% minVal = floor(min(ROIonly(:))/BinSize)*BinSize
% BinsCntrVect = [1, ((2:1:ceil(maxVal/BinSize))-0.5)*BinSize];
% ROIonlyUniform = round((BinSize-1)*(ROIonly-minVal)/(maxVal-minVal))+1;
% ROIonlyUniform = 1 + ceil((ceil((ROIonly-minVal)/BinSize) - 0.5) * BinSize);
% ROIonlyUniform2 = ceil(ROIonly/BinSize) - min(ceil(ROIonly/BinSize)) + 1; 

% The following is based on the definition of IBSI
% ROIonlyDiscr = ceil((ROIonly - minGL) / BinSize);
ROIonlyDiscr = floor((ROIonly - minGL) / BinSize)+1;
% ROIonlyUniform(ROIonlyUniform==0) = 1; 
ROIonlyUniform = ROIonlyDiscr; % - min(ROIonlyDiscr(:)) + 1; 
ROIonlyUniform(ROIonly==minGL) = 1; 
if min(ROIonlyUniform(:))<1, warning('Minimum discritization is below 1. Check the fixed bin size quantization.'); end 

levels = 1:max(ROIonlyUniform(:));
end