function SUV = getSUVmetrics(ROIonlyPET, pixelW,sliceS)
% -------------------------------------------------------------------------
% function [SUVmax,SUVpeak,SUVmean,aucCSH] = getSUVmetrics(ROIonlyPET)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes SUVmax, SUVpeak and SUVmean, AUC-CSH and Percent 
% Inactive metrics from the region of interest (ROI) of an input PET volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonlyPET: 3D array representing the PET volume in SUV format, with 
%               voxels outside the ROI set to NaNs. 
% -------------------------------------------------------------------------
% OUTPUTS:
% - SUVmax: Maximum SUV of the ROI.
% - SUVpeak: Average of the voxel with maximum SUV within the ROI and its 
%            26 connected neighbours.
% - SUVmean: Average SUV value of the ROI.
% - aucCSH: Area under the curve of the cumulative SUV-volume histogram
%           describing the percentage of total volume of the ROI above a 
%           percentage threshold of maximum SUV.
%           (van Velden et al., Eur J Nucl Med Mol Imaging 38(9), 2011).
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


% Initialization
ROIonlyPET = padarray(ROIonlyPET,[1 1 1],NaN);

% SUVmax
[SUV.SUVmax,indMax] = max(ROIonlyPET(:));

% SUVpeak (using 26 neighbors around SUVmax)
[indMaxX,indMaxY,indMaxZ] = ind2sub(size(ROIonlyPET),indMax);
connectivity = getneighbors(strel('arbitrary',conndef(3,'maximal')));
nPeak = length(connectivity);
neighborsMax = zeros(1,nPeak);
for i=1:nPeak
    neighborsMax(i) = ROIonlyPET(connectivity(i,1)+indMaxX,connectivity(i,2)+indMaxY,connectivity(i,3)+indMaxZ);
end
SUV.SUVpeak = mean(neighborsMax(~isnan(neighborsMax)));

% % % % SUVmax
% % % [SUVmax,indMax] = max(ROIonlyPET(:));
% % % 
% % % % SUVpeak (using 26 neighbors around SUVmax)
% % % % ################# TODO: THIS DOES NOT CHECK FOR MORE THAN ONE MAX VOXEL
% % % [indMaxX,indMaxY,indMaxZ] = ind2sub(size(ROIonlyPET),indMax);
% % % connectivity = getneighbors(strel('arbitrary',conndef(3,'maximal')));
% % % nPeak = length(connectivity);
% % % neighborsMax = zeros(1,nPeak);
% % % for i=1:nPeak
% % %     neighborsMax(i) = ROIonlyPET(connectivity(i,1)+indMaxX,connectivity(i,2)+indMaxY,connectivity(i,3)+indMaxZ);
% % % end
% % % SUVpeak = mean(neighborsMax(~isnan(neighborsMax)));

% % SUVpeak (using 26 neighbors around all voxeles in ROI,add by Wenbing Lv)
% [x,y,z] = size(ROIonlyPET);
% temp = 1;
% for i = 2:x-1
% for j = 2:y-1
% for k = 2:z-1
% cube26 = ROIonlyPET(i-1:i+1,j-1:j+1,k-1:k+1);
% meancube26 = mean(cube26(~isnan(cube26)));
% listall(temp) = meancube26;
% temp = temp+1;
% end
% end
% end
% SUVpeak = max(listall);


% SUVmean
SUV.SUVmean=mean(ROIonlyPET(~isnan(ROIonlyPET)));

%SUVstd
SUV.SUVstd = std(ROIonlyPET(~isnan(ROIonlyPET)));

% SUVvar 
SUV.SUVvar = (SUV.SUVstd)^2;

%SUVenergy
SUV.SUVenergy = sum(ROIonlyPET(~isnan(ROIonlyPET)).^2);

% AUC-CSH
SUV.AUC_CSH = getAUCCSH(ROIonlyPET);

% Skewness
SUV.Skewness = skewness(ROIonlyPET(:) , 1);

% Kurtosis
SUV.Kurtosis = kurtosis(ROIonlyPET(:) , 1);

% Median
SUV.Median = median(ROIonlyPET(:),'omitnan');

% SUVmin
SUV.SUVmin = min(ROIonlyPET(:),[],'omitnan');

% 10th percentile
SUV.Prcnt10 = prctile(ROIonlyPET(:),10);

% 90th percentile
SUV.Prcnt90 = prctile(ROIonlyPET(:),90);

% Interquantile range
SUV.IqntlRange = prctile(ROIonlyPET(:),75) - prctile(ROIonlyPET(:),25);

% Range
SUV.Range = SUV.SUVmax - SUV.SUVmin;

% Mean absolute deviation
nV = numel(find(~isnan(ROIonlyPET(:)))); % Number of voxels
SUV.MAD = sum(abs(ROIonlyPET(:) - SUV.SUVmean),'omitnan') / nV;

% Robust mean absolute deviation
RobustSet = find(ROIonlyPET(:)>SUV.Prcnt10 & ROIonlyPET(:)<SUV.Prcnt90);
nVr = nV - numel(RobustSet);
SUVmeanR = mean(ROIonlyPET(RobustSet));
SUV.RMAD = sum(abs(ROIonlyPET(RobustSet) - SUVmeanR),'omitnan') / nVr;

% Median absolute deviation
SUV.MedAD = sum(abs(ROIonlyPET(:) - SUV.Median),'omitnan') / nV;

% Coefficient of Variation
SUV.CoV = SUV.SUVstd / SUV.SUVmean;

% Quartile coefcient of dispersion
SUV.QCoD = (prctile(ROIonlyPET(:),75) - prctile(ROIonlyPET(:),25)) / (prctile(ROIonlyPET(:),75) + prctile(ROIonlyPET(:),25));

% Root Mean Square
SUV.RMS = sqrt(sum(ROIonlyPET(:).^2,'omitnan') / nV);








end