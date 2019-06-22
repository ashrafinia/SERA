function [ImgBoxResampQuntz3D,levels,MorphROI,IntsBoxROI,ImgWholeResmp,ROIwholeResmp,newpixelW,newsliceTh] = prepareVolume(volume,Mask,DataType,pixelW,sliceTh,R,newVoxelSize,VoxInterp,ROIInterp,ROI_PV,scaleType,isIsot2D,isScale,isGLrounding,DiscType,QuantAlg,Ng,isReSegRng,ResegIntrval,isOutliers)
% -------------------------------------------------------------------------
% function [ROIonly,levels,Maskout,ROIboxResmp,newpixelW,newsliceTh] =
% prepareVolume(volume, Mask, DataType, pixelW, sliceTh, R, newVoxelSize,
%               VoxInterp, ROIInterp, ROI_PV, scaleType, isIsot2D, isScale,
%               isGLrounding, DiscType, QuantAlg,Ng, isReSeg, ResegIntrval, 
%               isOutliers)    
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function prepares the input volume for 3D texture analysis. The 
% following operations are performed:
%
% 1. Computation of the smallest box containing region of interest (ROI), 
%    if necessary (ROIbox).
% 2. Pre-processing of the ROIbox (PET: square-root (default:off), MR:
%    Collewet normalizaton, CT: nothing).
% 3. Wavelet band-pass filtering (WBPF).
% 4. Isotropic resampling.
% 5. Image resegmentation.
% 6. grey levels rounding.
% 5. Quantization of intensity dynamic range.
%
% --> This function is compatible with both 2D and 3D analysis 
% -------------------------------------------------------------------------
% INPUTS:
% - volume: 2D or 3D array containing the medical images to analyze
% - mask: 2D or 3D array of dimensions corresponding to 'volume'. The mask 
%         contains 1's in the region of interest (ROI), and 0's elsewhere.
% - DataType: String specifying the type of scan analyzed. Either 'PET', 
%             'CT' or 'MRscan'.
% - pixelW: Numerical value specifying the in-plane resolution (mm) of 'volume'.
% - sliceTh: Numerical value specifying the slice spacing (mm) of 'volume'.
%           Put a random number for 2D analysis.
% - R: Numerical value specifying the ratio of weight to band-pass
%      coefficients over the weigth of the rest of coefficients (HHH and
%      LLL). Provide R=1 to not perform wavelet band-pass filtering.    
% - newVoxelSize: Numerical value specifying the scale at which 'volume' is
%                 isotropically  resampled in 3D (mm). 
% - VoxInterp: interpolation method for the intensity ROI
% - ROIInterp: interpolation method for the morphological ROI
% - ROI_PV: partial volume threshold for thresholding morphological ROI.
%           Used to threshold ROI after resampling: 
%           i.e. ROI(ROI<ROI_PV) =0, ROI(ROI>ROI_PV) = 1.  
% - scaleType: 'NoRescale' if no resampling, 'XYZscale' if resample in all
%              3 dimensions, 'XYscale' if only resampling in X and Y
%              dimensions, 'Zscale' if only scaling Z direction.
% - isIsot2D: =1 if resampling only in X and Y dimensions.
% - isScale: whether to perform resampling
% - isGLrounding: whether to perform grey level rounding
% - DiscType: discritization type. Either 'FNB' for fixed number of bins
%             or 'FBS' for fixed bin size. 
% - quantAlg: String specifying the quantization algorithm to use on
%             'volume'. Either 'Lloyd' for Lloyd-Max quantization, or
%             'Uniform' for uniform quantization. 
% - Ng: number of bins (for fixed number of bins method) or bin size
%            (for bin size method).
% - isReSeg: whether to perform resegmentation
% - ResegIntrval: a 1x2 array of values expressing resegmentation interval. 
% - isOutliers: whether to perform instensity outlier filtering 
% -------------------------------------------------------------------------
% OUTPUTS: ROIonly,levels,Maskout,ROIboxResmp,newpixelW,newsliceTh
% - ImgBoxResampQuntz3D: Smallest box containing the ROI, with the imaging
%           data of the ready for texture analysis computations. Voxels
%           outside the ROI are set to NaNs.  
% - levels: Vector containing the quantized gray-levels in the tumor region
%           (or reconstruction levels of quantization).
% - Maskout: Smallest matrix containing morphological ROI (0 or 1)
% - ROIboxResmp: Smallest matrix containing intensity ROI. This is the
%                resampled ROI, went through resegmtation and GL rounding,
%                and just before quantization. 
% - newpixelW: width of the voxel in the X (=Y) direction in resampled ROI
% - newsliceTh: Slice thickness of the voxel in the Z direction in resamped
%               ROI
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: July 2017
% - Revision Jan 2019: adding isclcGlobPeak to prevent resizing whole image
% and ROI only for the sake of global peak. 
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


%% VERIFICATION OF SOME INPUTS
% Set Quantization 
% minGL = min(volume(:)); 
if strcmp(DiscType,'FBS')
    quantization = @(x,y,z) fixedBinSizeQuantization(x,y,z);
elseif strcmp(DiscType,'FBN')
    quantization = @(x,y,z) uniformQuantization(x,y,z);
else
    error('Error with discretization type. Must either be "FBS" (Fixed Bin Size) or "FBN" (Fixed Number of Bins).')
end

if strcmp(QuantAlg,'Lloyd')
    quantization = @(x,y) lloydQuantization(x,y);
    %     elseif strcmp(quantAlgo,'Equal')
    %         quantization = @(x,y) equalQuantization(x,y);
    %     elseif strcmp(quantAlgo,'Uniform')
    %         quantization = @(x,y) uniformQuantization(x,y);
    %     elseif strcmp(quantAlgo,'FixedBinSize')
    %         quantization = @(x,y) fixedBinSizeQuantization(x,y);
    %     else
    %         error('Error with quantization algorithm input. Must either be ''Equal'' or ''Lloyd'' or ''Uniform''')
end



%% COMPUTATION OF THE SMALLEST BOX CONTAINING THE ROI
% [boxBound] = computeBoundingBox(Mask);
% maskBox = Mask(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
% ROIbox = volume(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
ROIBox = Mask;
Imgbox = volume;

%% PRE-PROCESSING OF ROI BOX
Imgbox = double(Imgbox);
if strcmp(DataType,'MRscan')
    ROIonly = Imgbox;
    ROIonly(~ROIBox) = NaN;
    temp = CollewetNorm(ROIonly);
    ROIBox(isnan(temp)) = 0;
% elseif strcmp(DataType,'PETscan')
%    ROIbox = sqrt(ROIbox);
end


%% WAVELET BAND-PASS FILTERING
if R ~= 1
    if sum(isnan(Imgbox(:)))
        Imgbox = fillBox(Imgbox); % Necessary in cases we have a ROI box containing NaN's.
    end
    Imgbox = waveletBPfilt(Imgbox,R,'sym8');
end


%% ISOTROPIC RESAMPLING
flagPW = 0;
if strcmp(scaleType,'NoRescale')
    flagPW = 0;
elseif strcmp(scaleType,'XYZscale')
    flagPW = 1;
elseif strcmp(scaleType,'XYscale')
    flagPW = 2;
elseif strcmp(scaleType,'Zscale')
    flagPW = 3;
end

if isIsot2D
    flagPW = 2;
end

if isScale == 0
    flagPW = 0; 
end

if flagPW == 0
    a = 1;
    b = 1; 
    c = 1;
elseif flagPW ==1
    a = pixelW/newVoxelSize;
    b = pixelW/newVoxelSize;
    c = sliceTh/newVoxelSize;
elseif flagPW == 2
    a = pixelW/newVoxelSize;
    b = pixelW/newVoxelSize;
    c = 1;
elseif flagPW == 3
    a = 1;
    b = 1;
    c = sliceTh/pixelW;
end

%% Resampling
ImgBoxResmp = Imgbox; ImgWholeResmp = volume;
ROIBoxResmp = ROIBox; ROIwholeResmp = Mask;
if numel(size(Imgbox))==3 && (flagPW~=0)
    if a + b + c ~= 3 % If false, no resampling is needed
        % intensity mask and morphological mask
        ROIBoxResmp = imresize3D(ROIBox,[a, b, c],[ceil(size(ROIBox,1)*a),ceil(size(ROIBox,2)*b),ceil(size(ROIBox,3)*c)],ROIInterp ,'fill');
        Imgbox(isnan(Imgbox)) = 0; 
        ImgBoxResmp = imresize3D(Imgbox,[a, b, c],[ceil(size(Imgbox,1)*a),ceil(size(Imgbox,2)*b),ceil(size(Imgbox,3)*c)],VoxInterp,'fill');
        ROIBoxResmp(ROIBoxResmp<ROI_PV) = 0; 
        ROIBoxResmp(ROIBoxResmp>=ROI_PV) = 1;  
        % Raw image and mask for SUVpeak
        ROIwholeResmp = imresize3D(Mask,[a, b, c],[ceil(size(Mask,1)*a),ceil(size(Mask,2)*b),ceil(size(Mask,3)*c)],ROIInterp,'fill');
        ImgWholeResmp = imresize3D(volume,[a, b, c],[ceil(size(volume,1)*a),ceil(size(volume,2)*b),ceil(size(volume,3)*c)],VoxInterp ,'fill');
        if max(ROIwholeResmp(:))<ROI_PV
            disp('Resampled ROI has no voxels with value above ROI_PV. Cutting ROI_PV to half.');
            ROI_PV = ROI_PV/2;
        end
        ROIwholeResmp(ROIwholeResmp<ROI_PV) = 0;
        ROIwholeResmp(ROIwholeResmp>=ROI_PV) = 1;
    end
elseif numel(size(Imgbox))==2 && (flagPW~=0)
    if a + b ~= 2 % If false, no resampling is needed
        % intensity mask and morphological mask
        ROIBoxResmp = imresize(ROIBox,[ceil(size(ROIBox,1)*a),ceil(size(ROIBox,2)*b)],ROIInterp);
        ImgBoxResmp = imresize(Imgbox,[ceil(size(Imgbox,1)*a),ceil(size(Imgbox,2)*b)],VoxInterp ,'Antialiasing',true);
        ROIBoxResmp(ROIBoxResmp<ROI_PV) = 0; 
        ROIBoxResmp(ROIBoxResmp>=ROI_PV) = 1;
        %         ROIbox = ROIbox .* maskBox;
        % Raw image and mask for SUVpeak
        ROIwholeResmp = imresize(Mask,[ceil(size(Mask,1)*a),ceil(size(Mask,2)*b)],ROIInterp);
        ImgWholeResmp = imresize(volume,[ceil(size(volume,1)*a),ceil(size(volume,2)*b)],VoxInterp ,'Antialiasing',true);
        if max(ROIwholeResmp(:))<ROI_PV
            disp('Resampled ROI has no voxels with value above ROI_PV. Cutting ROI_PV to half.');
            ROI_PV = ROI_PV/2;
        end
        ROIwholeResmp(ROIwholeResmp<ROI_PV) = 0;
        ROIwholeResmp(ROIwholeResmp>=ROI_PV) = 1;
    end
end

ImgBoxResmp(~ROIBoxResmp) = NaN;

%% GL rounding 
% IntsBoxROI = roundGL(ImgBoxResmp , isGLrounding); 
% ImgWholeResmp = roundGL(ImgWholeResmp , isGLrounding); 

IntsBoxROI = roundGL(ImgBoxResmp , isGLrounding); 
ImgWholeResmp = roundGL(ImgWholeResmp , isGLrounding); 

IntsBoxROItmp1 = IntsBoxROI; 
ImgWholeResmptmp1 = ImgWholeResmp; 
IntsBoxROItmp2 = IntsBoxROI; 
ImgWholeResmptmp2 = ImgWholeResmp; 
%% Image resegmentation
% 1- Range re-segmentation
if isReSegRng
    IntsBoxROItmp1(IntsBoxROI<ResegIntrval(1)) = NaN;
    IntsBoxROItmp1(IntsBoxROI>ResegIntrval(2)) = NaN;
    ImgWholeResmptmp1(ImgWholeResmp<ResegIntrval(1)) = NaN;
    ImgWholeResmptmp1(ImgWholeResmp>ResegIntrval(2)) = NaN;
end

% Intensity outlier filtering
if isOutliers
    Mu = mean(IntsBoxROI(:),'omitnan');
    Sigma = std(IntsBoxROI(:),'omitnan');
    IntsBoxROItmp2(IntsBoxROI<(Mu-3*Sigma)) = NaN;
    IntsBoxROItmp2(IntsBoxROI>(Mu+3*Sigma)) = NaN;
        % Raw image and mask for SUVpeak
    Mu = mean(ImgWholeResmp(:),'omitnan');
    Sigma = std(ImgWholeResmp(:),'omitnan');
    ImgWholeResmptmp2(ImgWholeResmp<(Mu-3*Sigma)) = NaN;
    ImgWholeResmptmp2(ImgWholeResmp>(Mu+3*Sigma)) = NaN;
end

% Get the intersection of both re-segmentations 
IntsBoxROI      = getMutualROI(IntsBoxROItmp1 , IntsBoxROItmp2);
ImgWholeResmp   = getMutualROI(ImgWholeResmptmp1 , ImgWholeResmptmp2); 

%% New slice thickness and pixel width
newpixelW = pixelW / a;
newsliceTh = sliceTh / c; 


%% Determing minimum GL for discretization
if strcmp(DataType, 'PET')
    minGL = 0; 
elseif strcmp(DataType, 'CT')
    if isReSegRng
%         minGL = max(min(IntsROI(:)), ResegIntrval(1));  
        minGL = ResegIntrval(1);
    else
        minGL = min(IntsBoxROI(:)); 
    end
else
    minGL = min(IntsBoxROI(:)); 
end


%% QUANTIZATION
% Run quantization
[ImgBoxResampQuntz3D,levels] = quantization(IntsBoxROI,Ng,minGL);


%% Crop the ROI to get rid of extra zeros
[boxBound] = computeBoundingBox(ROIBoxResmp);
MorphROI = ROIBoxResmp(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
% [boxBound] = computeBoundingBox(IntsROI); %%%% THIS NEEDS TO BE VERIFIED, IF WE NEED TO HAVE THIS LINE OR NOT.#########
IntsBoxROI = IntsBoxROI(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
ImgBoxResampQuntz3D = ImgBoxResampQuntz3D(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

% if boxBound(1,2)==boxBound(1,1)
%     MorphROI = padarray(MorphROI,[1 1 1]);
%     IntsBoxROI = padarray(IntsBoxROI,[1 1 1]);
%     ImgBoxResampQuntz3D = padarray(ImgBoxResampQuntz3D,[1 1 1]);
% end

end






