function [AllFeats] = Calc_All_Features(RawImg,ROI,pixelW,sliceTh,IsotVoxSize,IsotVoxSize2d,qntzAlg,bin,DiscType,DataType,VoxInterp,ROIInterp,ROI_PV,isIsot2D,isScale,isGLrounding,isReSeg,ResegIntrval,isQuntzStat,IVHconfig,isOutliers,Feats2out)
% -------------------------------------------------------------------------
% [AllFeats] = Calc_All_Features(RawImg, ROI, pixelW, sliceTh, IsotVoxSize,
%                                IsotVoxSize2d, qntzAlg, bin, DiscType,
%                                DataType, VoxInterp, ROIInterp, ROI_PV,
%                                isIsot2D, isScale, isGLrounding, isReSeg,
%                                ResegIntrval, isQuntzStat, isOutliers)    
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function prepares the input volume for 2D and 3D texture analysis.
% Then calculates various radiomics features, including: morphological,
% first-order: statistical, histogram, volume histogram; second-order: GLCM
% and GLRLM; and higher order: GLSZM, GLDZM, NGLDM, NGTDM.
% -------------------------------------------------------------------------
% INPUTS:
% - RawImg: the 2D or 3D matrix of intensities
% - ROI: the 2D or 3D matrix of the mask. Has to be the same size as RawImg
% - pixelW: Numerical value specifying the in-plane resolution(mm) of RawImg
% - sliceTh: Numerical value specifying the slice spacing (mm) of RawImg
%           Put a random number for 2D analysis.
% - IsotVoxSize: Numerical value specifying the new voxel size (mm) at
%                which the RIO is isotropically  resampled to in 3D.
% - IsotVoxSize2D: Numerical value specifying the new voxel size (mm) at
%                which the RIO is isotropically  resampled to in 3D. 
% - quantAlg: String specifying the quantization algorithm to use on
%             'volume'. Either 'Lloyd' for Lloyd-Max quantization, or
%             'Uniform' for uniform quantization. 
% - bin: number of bins (for fixed number of bins method) or bin size
%        (for bin size method). It can be an array of bins, and the code
%        will loop over each bin. 
% - DiscType: discritization type. Either 'FNB' for fixed number of bins
%             or 'FBS' for fixed bin size. 
% - DataType: String specifying the type of scan analyzed. Either 'PET', 
%             'CT', or 'MRscan'.
% - VoxInterp: interpolation method for the intensity ROI
% - ROIInterp: interpolation method for the morphological ROI
% - ROI_PV: partial volume threshold for thresholding morphological ROI.
%           Used to threshold ROI after resampling: 
%           i.e. ROI(ROI<ROI_PV) =0, ROI(ROI>ROI_PV) = 1.  
% - isIsot2D: =1 if resampling only in X and Y dimensions.
% - isScale: whether to perform resampling
% - isGLrounding: whether to perform grey level rounding
% - isReSeg: whether to perform resegmentation
% - ResegIntrval: a 1x2 array of values expressing resegmentation interval. 
% - isQuntzStat: whether to use quantized image to calculate first order
%                statistical features. If 0, no image resample/interp for
%                calculating statistical features. (0 is preferrable for
%                PET images).
% - isOutliers: whether to perform instensity outlier filtering.
% - Feats2out: the option of which calculated radiomic features to return.
% -------------------------------------------------------------------------
% OUTPUTS: 
% - AllFeats: A vector (for a single bin) or a matrix (for multiple bins)
%             of calculated Radiomic features.
% -------------------------------------------------------------------------
% % AUTHOR(S): 
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2017
% - Revision: Nov 2018
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of Radiomics Package by Saeed Ashrafinia, Rahmimlab.com
% --> Copyright (C) 2013-2017  Saeed Ashrafinia, Johns Hopkins University
% 
%   All rights reserved for Saeed Ashrafinia and Arman Rahmim. 
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% -------------------------------------------------------------------------

%% Initialization
Lbin = numel(bin);
img = RawImg;
img(find(ROI==0)) = NaN; %#ok<FNDSB>

%% TextureFeatures for different bin number
AllFeats = []; 
for m = 1:Lbin
    %% ROIs preperation 
    % 1. General 3D ROI:
    % This intensity ROI is for first order and Morphological features.
    [ROIBox3D,levels3D,ROIonlyMorph3D,IntsROI,RawInts,RawROI,newPixW,newSliceTh] = prepareVolume(RawImg,ROI,DataType,pixelW,sliceTh,1,IsotVoxSize,VoxInterp,ROIInterp,ROI_PV,'XYZscale',isIsot2D,isScale,isGLrounding,DiscType,qntzAlg,bin(m),isReSeg,ResegIntrval,isOutliers);
    [ImgBox, ~] = getImgBox(img,ROI,isReSeg,ResegIntrval);
    
    % 2. 2D and 3D resampled ROI for higher order features
    % If "isotropic 2D" is selected, then ROI_2D = ROI_3D = ROIBoxG above.
    if isIsot2D || pixelW == IsotVoxSize || ~isScale
        ROIBox2D = ROIBox3D;  levels2D = levels3D; ROIonly2D = ROIonlyMorph3D;
    else
        [ROIBox2D,levels2D,ROIonly2D] = prepareVolume(RawImg,ROI,DataType,pixelW,sliceTh,1,IsotVoxSize2d,VoxInterp,ROIInterp,ROI_PV,'XYscale',0,isScale,isGLrounding,DiscType,qntzAlg,bin(m),isReSeg,ResegIntrval,isOutliers);
    end
    
    %% (1) Morphological features
    [MorphVect] = getMorph(ROIBox3D,ROIonlyMorph3D,IntsROI, newPixW,newSliceTh);
    
    %% (1.5) SUVpeak (local intensity) calculation
    [SUVpeak] = getSUVpeak(RawInts,RawROI,newPixW,newSliceTh);
    
    %% (2) Statistical features
    % whether to use the above quantized image to calculate statistical
    % features, or just use the raw image and original voxel size. 
    % For PET images, set isQuntzStat=0.  
    if isQuntzStat      
        [StatsVect] = getStats(IntsROI);
    else
        [StatsVect] = getStats(ImgBox);
    end
    
    %% (3) Moment Invariants
    MI_feats = getMI(ImgBox)';
    
    %% (4) and (5) Histogram and Intensity Histogram features
    [HistVect]  = getHist(ROIBox3D,bin(m), DiscType);
    [IVHvect]   = getIntVolHist(IntsROI,ROIBox3D,bin(m),isReSeg,ResegIntrval,IVHconfig);
    
    %% (6) GLCM
    [GLCM2D_KSKD, GLCM2D_MSKD, GLCM2D_KSMD, GLCM2D_MSMD] = getGLCM2Dtex(ROIBox2D,levels2D); 
    [GLCM3D_Cmb, GLCM3D_Avg] = getGLCM3Dtex(ROIBox3D,levels3D);
    
    %% (7) GLRLM
    [GLRLM2D_KSKD, GLRLM2D_KSMD, GLRLM2D_MSKD, GLRLM2D_MSMD] = getGLRLM2Dtex(ROIBox2D,levels2D);
    [GLRLM3D_Cmb, GLRLM3D_Avg] = getGLRLM3Dtex(ROIBox3D,levels3D);
    
    %% (8) GLSZM
    [GLSZM2D, GLSZM3D, GLSZM25D ] = getGLSZMtex(ROIBox2D,ROIBox3D,levels2D,levels3D);
    
    %% (9) NGTDM
    [NGTDM2D, NGTDM3D, NGTDM25D] = getNGTDMtex(ROIBox2D,ROIBox3D,levels2D,levels3D);
    
    %% (10) NGLDM
    [NGLDM2D, NGLDM3D, NGLDM25D] = getNGLDMtex(ROIBox2D,ROIBox3D,levels2D,levels3D);
    
    %% (11) GLDZM
    [GLDZM2D, GLDZM3D, GLDZM25D] = getGLDZMtex(ROIBox2D,ROIBox3D,ROIonly2D,ROIonlyMorph3D,levels2D,levels3D);
    
    %% Save all
    [AllFeats] = ReturnFeatures(AllFeats, Feats2out , MorphVect, SUVpeak, StatsVect, HistVect, IVHvect, ...
        GLCM2D_KSKD,  GLCM2D_KSMD,  GLCM2D_MSKD,  GLCM2D_MSMD,  GLCM3D_Avg,  GLCM3D_Cmb, ...
        GLRLM2D_KSKD, GLRLM2D_KSMD, GLRLM2D_MSKD, GLRLM2D_MSMD, GLRLM3D_Avg, GLRLM3D_Cmb, ...
        GLSZM2D, GLSZM25D, GLSZM3D, GLDZM2D, GLDZM25D, GLDZM3D, ...
        NGTDM2D, NGTDM25D, NGTDM3D, NGLDM2D, NGLDM25D, NGLDM3D, ...
        MI_feats);
    
end


end





