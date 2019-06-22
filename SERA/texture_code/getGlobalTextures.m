function [HistVect] = getGlobalTextures(ROIonly,BinSize, DiscType,DataType,minGL)
% -------------------------------------------------------------------------
% function [textures] = getGlobalTextures(ROIonly,Nbins)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes Global texture features from the region of 
% interest (ROI) of an input volume.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: Smallest box containing the ROI, with the imaging data ready 
%            for texture analysis computations. Voxels outside the ROI are 
%            set to NaNs.
% - levels: Vector containing the quantized gray-levels in the tumor region
%           (or reconstruction levels of quantization).
%
% ** 'ROIonly' and 'levels' should be outputs from 'prepareVolume.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different Global texture
%             features as defined below.
% -------------------------------------------------------------------------
% AUTHOR(S):    Saeed Ashrafinia
%               Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013.
% - Revision: July 2017
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2013-2015  Martin Vallieres
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


% PRELIMINARY
% vectorValid = ROIonly(~isnan(ROIonly)); % Changed this to the below to avoid an error caused by calculating "u" below for a single-dimension-Z vector
ROIarrayValid = squeeze(ROIonly(~isnan(ROIonly)));

if strcmp(DiscType, 'FNB')% strcmp(DataType, 'PET') &&
    % Fixed number of bins
    histraw = hist(ROIarrayValid,BinSize);
    histo = histraw./(sum(histraw(:)));
    vectNg = 1:BinSize;
    % try
    Mean = histo*vectNg';
    % catch
    %     u = 0;
    %     warning('Small ROI; set u=0 in getGlobalTexture.m');
    % end
elseif strcmp(DiscType, 'FBS') %strcmp(DataType, 'CT') || strcmp(DataType, 'MRI')
    % Fixed bin size
%     BinsCntrVect = [1, ((2:1:ceil(max(vectorValid(:))/BinSize))-0.5)]*BinSize; % 1:length(vectorValid); %
    BinsCntrVect = ceil(minGL:BinSize:max(ROIarrayValid(:)));
    vectNg = BinsCntrVect;
%     histraw = hist(ROIarrayValid,BinsCntrVect);
    [histraw]=histcounts(ROIarrayValid,[BinsCntrVect, (BinsCntrVect(end)+BinSize)]);
    histo = histraw./(sum(histraw(:)));
    BinSize = numel(vectNg);
%     Mean = histo*vectNg';
    Binslist = 1:BinSize;
    Mean = histo * Binslist';
end





% COMPUTATION OF TEXTURES
% 1. Mean
Histogram.Mean = Mean;

% 2. Variance
Var = 0;
for i=1:BinSize
    Var = Var+histo(i)*(i-Mean)^2;
end
sigma = sqrt(Var);
Histogram.Variance = Var;

% 3. Skewness
skw = 0;
for i = 1:BinSize
    skw = skw+histo(i)*(i-Mean)^3;
end
skw = skw/sigma^3;
Histogram.Skewness = skw;
% textures.Skewness = skewness(vectorValid,1);

% 4. Kurtosis
krt = 0;
for i = 1:BinSize
    krt = krt+histo(i)*(i-Mean)^4;
end
krt = (krt/sigma^4) - 3;
Histogram.Kurtosis = krt;
% textures.Kurtosis = kurtosis(vectorValid,1);

% Median
Histogram.Median = find(cumsum(histo)>=0.5,1,'first');

% Hist Min
Histogram.SUVmin = min(Binslist);

% Hist Max
Histogram.SUVmax = max(Binslist);

% 10th percentile
Histogram.Prcnt10 = find(cumsum(histo)>=0.1,1,'first');

% 90th percentile
Histogram.Prcnt90 = find(cumsum(histo)>=0.9,1,'first');

% Intensity Histogram Mode: The mode of Xd is the most common discretised grey level present
[~, ModeIdx] = max(histo); 
[~,Idx] = sort(abs(ModeIdx-Mean)); 
Histogram.Mode = ModeIdx(Idx);

% Interquantile range
Histogram.IqntlRange = find(cumsum(histo)>=0.75,1,'first') - find(cumsum(histo)>=0.25,1,'first');

% Range
Histogram.Range = Histogram.SUVmax - Histogram.SUVmin;

% Mean absolute deviation
nV = numel(find(~isnan(ROIarrayValid))); % Number of voxels
Histogram.MAD = sum(abs(histraw.*(Binslist - Mean)))/nV;

% Robust mean absolute deviation
RobustSet = find(Binslist>=Histogram.Prcnt10 & Binslist<=Histogram.Prcnt90);
nVr = sum(histraw(RobustSet));
SUVmeanR = histraw(RobustSet)./nVr * Binslist(RobustSet)'; %mean(Binslist(RobustSet));
% Histogram.RMAD = sum(abs(Binslist(RobustSet) - SUVmeanR)) / nVr;
Histogram.RMAD = sum(abs(histraw(RobustSet).*(Binslist(RobustSet) - SUVmeanR))) / nVr;

% Median absolute deviation
% Histogram.MedAD = sum(abs(ROIarrayValid - Histogram.Median),'omitnan') / nV;
Histogram.MedAD = sum(abs(histraw.*(Binslist - Histogram.Median)))/nV;

% Coefficient of Variation
Histogram.CoV = sqrt(Var) / Mean;

% Quartile coefcient of dispersion
Histogram.QCoD = (find(cumsum(histo)>=0.75,1,'first') - find(cumsum(histo)>=0.25,1,'first')) / (find(cumsum(histo)>=0.75,1,'first') + find(cumsum(histo)>=0.25,1,'first'));


% Energy (aka Uniformity)
Histogram.Energy = sum(histo.^2);

% Entropy
entropy = 0;
for i=1:BinSize
    entropy = entropy-histo(i)*log2(histo(i)+realmin);
end
Histogram.Entropy = entropy;

% Histogram Grandient based features:
% H = [];
% Maximum histogram gradient


HistVect = [Histogram.Mean; Histogram.Variance; Histogram.Skewness; Histogram.Kurtosis; Histogram.Median;...
            Histogram.SUVmin; Histogram.Prcnt10; Histogram.Prcnt90; Histogram.SUVmax; Histogram.Mode; ...
            Histogram.IqntlRange; Histogram.Range; Histogram.MAD; Histogram.RMAD; Histogram.MedAD; ...
%             Histogram.CoV; Histogram.QCoD; Histogram.Entropy; Histogram.Energy;  [NaN; NaN; NaN; NaN];  ];
            Histogram.CoV; Histogram.QCoD; Histogram.Entropy; Histogram.Energy;    ];



end






