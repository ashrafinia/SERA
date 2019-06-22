function [ROIonlyUniform,levels] = uniformQuantization(ROIonly,Ng,minGL)
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
% - Revision: July 2017
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
maxGL = max(ROIonly(:));
% minGL = min(ROIonly(:));
% ROIonly(ROIonly<minGL) = minGL; 
% ROIonlyUniform = round((Ng-1)*(ROIonly-minGL)/(maxGL-minGL))+1;

ROIonlyUniform = ceil((Ng)*(ROIonly-minGL)/(maxGL-minGL));
ROIonlyUniform(ROIonlyUniform==0)=1;

if min(ROIonlyUniform(:)) < 1 
    error('Something wrong with FNB discritization.'); 
end

levels = 1:Ng;
end