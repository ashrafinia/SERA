function [MorphVect] = getMorph(MorphBox,MorphROI,ROIints, pixelW,sliceS)
% -------------------------------------------------------------------------
% function [MorphVect] = getMorph(MorphBox,MorphROI,ROIints, pixelW,sliceS)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates morphological features.
% -------------------------------------------------------------------------
% INPUTS:
% - MorphBox: The smallest box containing the resampled 3D ROI, with the
%             imaging data ready for texture analysis computations. Voxels
%             outside the ROI are set to NaNs.
% - MorphROI: The smallest box containing the 3D morphological ROI. Voxels
%             outside the ROI are set to 0. It has only 0 and 1.
% - ROIints: The smallest box containing the 3D intensity ROI with
%             their intensities. Voxels outside the ROI are set to NaNs.
% - pixelW: width of the voxel in the X (=Y) direction
% - sliceS: Slice thickness of the voxel in the Z direction
% 
% Note: The first 3 parameters are outputs of prepareVolume.m
% -------------------------------------------------------------------------
% OUTPUTS:
% - MorphVect: An array of calculated morphological features
% -------------------------------------------------------------------------
% AUTHOR(S): Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2017
% - Modification: July 2017
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
% Set values of img outside of ROI to NaN; necessary for later processing
% (e.g. in getSUVmetrics, when mean is calculated only for non-NaN voxels)
% % MorphROI = zeros(size(MorphBox));
% % MorphROI(~isnan(MorphBox)) = 1;
% MorphROI(MorphROI==0) = NaN;

[ROIboxF,ROIboxV] = getMesh(MorphROI,pixelW,pixelW,sliceS); % Convert ROI to mesh
[MV, Surface] = getVolSurfaceMesh(ROIboxF,ROIboxV);

MorphROI(MorphROI==0) = NaN;
[MVappx] = getVolume(MorphROI,pixelW,sliceS);% /1e3; we suppressed it for the current project

[xM,yM,zM]=ind2sub(size(MorphROI),find(~isnan(MorphROI))); 
[xI,yI,zI]=ind2sub(size(ROIints),find(~isnan(ROIints))); 
try
CoM_shift = norm(mean([xM*pixelW yM*pixelW zM*sliceS],1,'omitnan') - ...
            sum([xI*pixelW yI*pixelW zI*sliceS].*repmat(ROIints(find(~isnan(ROIints))),[1,3]),1)./repmat(sum(ROIints(:),1,'omitnan'),[1,3])); %#ok<FNDSB>
catch
    disp('**** ROI is 1D. Brace for problems.....');
    CoM_shift = abs(mean(squeeze(MorphROI)-squeeze(ROIints))); 
end

%MarginSharpness = getMarginSharpness(ROI,IntSUVPT);
%% Compute the shape features
% Shape = getShape(MV*1e3,Surface);
[Shape] = getShape(MV,Surface);
compactness1 = Shape.compactness1;  % Compactness 1 is a measure for how compact, or sphere-like the volume is.
compactness2 = Shape.compactness2;  % Compactness 2 is another measure to describe how sphere-like the volume is.
sphericity = Shape.sphericity;      % Sphericity is a further measure to describe how sphere-like the volume is.
SVratio = Shape.SVratio;
%Irregularity = Shape.Irregularity;
sphericalDisp = Shape.sphericalDisp;% Spherical disproportion is another measure to describe how sphere-like the volume is.
Asphericity = Shape.Asphericity;    % Asphericity describes how much the ROI deviates from a perfect sphere.

% [Eccentricity] = getEccentricity(MorphROI,pixelW,sliceS);
% thresh = 0.005*max(max(max(MorphBox)));
% [PercentInactive] = getPercentInactive(MorphROI,thresh);
% [SizeROI] = getSize(SUVbox,pixelW,sliceS);

% [~ , CHarea] = getSolidity(SUVbox,pixelW,sliceS); % Solidity is also Volume Density of Convex Hull
% AreaDensConvHull = Surface / CHarea;


%% Volume Densities 
[VDs, ConvHullV] = getVDs(MorphROI,ROIboxV,pixelW,sliceS,MV,Surface );
Mean = mean(ROIints(:),'omitnan');
TLG = Mean * MV;
Max3Ddiam = getMax3Ddiam(ROIboxV, ConvHullV);

%% Autocorrelation metrics
MoranI = NaN; GearyC = NaN; 
[MoranI, GearyC] = getAutoCorrs(ROIints, xI,yI,zI , pixelW,sliceS , Mean); % THIS CODE IS NOT WORKING. 

%% GT2: Defined by Nikolay from UBC; sums up volume of uptake with uptake > 2 (threshold)
threshold=2;
% [GT2] = getVolume_gt(MorphBox,threshold,pixelW,sliceS);

%% Export
MorphVect = [MV; MVappx; Surface; SVratio; compactness1; compactness2; ...
            sphericalDisp; sphericity; Asphericity; CoM_shift; Max3Ddiam; ...
            VDs(1:end); TLG; MoranI; GearyC ];  %; ...
%             Eccentricity; PercentInactive; GT2];  % The last three are not in the standardization project
        
        
end
% VDs(1:7); [NaN; NaN]; VDs(8:13); TLG; [NaN; NaN]; ...
        