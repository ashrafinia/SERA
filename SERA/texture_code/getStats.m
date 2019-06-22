function [StatsVect] = getStats(ROIbox)
% -------------------------------------------------------------------------
% function [SUVmax,SUVpeak,SUVmean,aucCSH] = getSUVmetrics(ROIonlyPET)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes SUVmax, SUVpeak and SUVmean, AUC-CSH and Percent 
% Inactive metrics from the region of interest (ROI) of an input PET volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIbox: The smallest box containing the resampled 3D ROI, with the
%           imaging data ready for texture analysis computations. Voxels
%           outside the ROI are set to NaNs.
% -------------------------------------------------------------------------
% OUTPUTS:
% A list of 18 statistical features as documented by ISBI. 
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


% Initialization
ROIboxPadded = padarray(ROIbox,[1 1 1],NaN);

% SUVmax
StatsStr.SUVmax = max(ROIboxPadded(:));

% SUVmean
StatsStr.SUVmean=mean(ROIboxPadded(~isnan(ROIboxPadded)));

%SUVstd
StatsStr.SUVstd = std(ROIboxPadded(~isnan(ROIboxPadded)));

% SUVvar 
StatsStr.SUVvar = var(ROIboxPadded(~isnan(ROIboxPadded)),1);

%SUVenergy
StatsStr.Energy = sum(ROIboxPadded(~isnan(ROIboxPadded)).^2);

% Skewness
StatsStr.Skewness = skewness(ROIboxPadded(:) , 1);

% Kurtosis
StatsStr.Kurtosis = kurtosis(ROIboxPadded(:) , 1) - 3;

% Median
StatsStr.Median = median(ROIboxPadded(:),'omitnan');

% SUVmin
StatsStr.SUVmin = min(ROIboxPadded(:),[],'omitnan');

% 10th percentile
StatsStr.Prcnt10 = prctile(ROIboxPadded(:),10);

% 90th percentile
StatsStr.Prcnt90 = prctile(ROIboxPadded(:),90);

% Interquantile range
StatsStr.IqntlRange = prctile(ROIboxPadded(:),75) - prctile(ROIboxPadded(:),25);

% Range
StatsStr.Range = StatsStr.SUVmax - StatsStr.SUVmin;

% Mean absolute deviation
nV = numel(find(~isnan(ROIboxPadded(:)))); % Number of voxels
StatsStr.MAD = sum(abs(ROIboxPadded(:) - StatsStr.SUVmean),'omitnan') / nV;

% Robust mean absolute deviation
RobustSet = find(ROIboxPadded(:)>=StatsStr.Prcnt10 & ROIboxPadded(:)<=StatsStr.Prcnt90);
nVr = numel(RobustSet);
SUVmeanR = sum(ROIboxPadded(RobustSet))/ nVr;
StatsStr.RMAD = sum(abs(ROIboxPadded(RobustSet) - SUVmeanR),'omitnan') / nVr;

% Median absolute deviation
StatsStr.MedAD = sum(abs(ROIboxPadded(:) - StatsStr.Median),'omitnan') / nV;

% Coefficient of Variation
StatsStr.CoV = sqrt(StatsStr.SUVvar) / StatsStr.SUVmean;

% Quartile coefcient of dispersion
StatsStr.QCoD = (prctile(ROIboxPadded(:),75) - prctile(ROIboxPadded(:),25)) / (prctile(ROIboxPadded(:),75) + prctile(ROIboxPadded(:),25));

% Root Mean Square
StatsStr.RMS = sqrt(sum(ROIboxPadded(:).^2,'omitnan') / nV);

% % AUC-CSH
% StatsStr.AUC_CSH = getAUCCSH(SUVboxPadded);



% export
StatsVect = [StatsStr.SUVmean;  StatsStr.SUVvar;     StatsStr.Skewness;   StatsStr.Kurtosis; ...
            StatsStr.Median;   StatsStr.SUVmin;     StatsStr.Prcnt10;    StatsStr.Prcnt90;...
            StatsStr.SUVmax;   StatsStr.IqntlRange; StatsStr.Range;      StatsStr.MAD; ...
            StatsStr.RMAD;     StatsStr.MedAD;      StatsStr.CoV;        StatsStr.QCoD; ...
            StatsStr.Energy;   StatsStr.RMS;        ];

        


end