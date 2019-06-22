function [NGLDM] = getNGLDM(ROIonly,levels)
% -------------------------------------------------------------------------
% function [NGLDM] = getNGLDM(ROIonly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the Neighbouring grey level dependence Matrix
% (NGLDM) of the region of interest (ROI) of an input volume. The input
% volume is assumed to be isotropically resampled. Only one NGLDM is
% computed per scan.    
%
% --> This function is compatible with 2D analysis 
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Sun, C. and Wee, W. G. (1983). Neighboring gray level dependence
%       matrix for texture classification. 
% [2] Image biomarker standardisation initiative work document
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
% - NGLDM: Gray-Level Run-Length Matrix of 'ROIonly'.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Saeed Ashrafinia 
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2017
% -------------------------------------------------------------------------
% --> Copyright (c) 2007-2017, Saeed Ashrafinia
%     All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------

%INITIALIZATION
a = 0;  % Positive integer coarseness parameter
d = 1;  % distance of the neighborhood


% PRELIMINARY
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
% QUANTIZATION EFFECTS CORRECTION
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
% uniqueVol = round(levels*adjust)/adjust;
ROIonly = round(ROIonly*adjust)/adjust;
sizeV = size(ROIonly);


% START COMPUTATION
if numel(sizeV) == 3
    nComp = sizeV(3); 
    NnInit = (2*d+1)^3;
    ROIonly = padarray(ROIonly,[1,1,1],NaN,'both');
else
    nComp = 1; 
    NnInit = (2*d+1)^2;
    ROIonly = padarray(ROIonly,[1,1],NaN,'both');
end
NGLDM = zeros(nLevel,NnInit);



% COMPUTATION OF NGTDM
if nComp == 1
    [I,J] = ind2sub(size(ROIonly),find(~isnan(ROIonly)));
    for n = 1:length(I)
        XglC = ROIonly(I(n),J(n));
        k = numel(find(abs(XglC - ROIonly((I(n)-1):(I(n)+1) , (J(n)-1):(J(n)+1))) <= a));
        NGLDM(XglC,k) = NGLDM(XglC,k) + 1;
    end
else
    [I,J,K] = ind2sub(size(ROIonly),find(~isnan(ROIonly)));
    for n = 1:length(I)
        XglC = ROIonly(I(n),J(n),K(n));
        k = numel(find(abs(XglC - ROIonly((I(n)-1):(I(n)+1) , (J(n)-1):(J(n)+1) , (K(n)-1):(K(n)+1))) <= a));
        NGLDM(XglC,k) = NGLDM(XglC,k) + 1;
    end
end


% REMOVE UNECESSARY COLUMNS
stop = find(sum(NGLDM),1,'last');
NGLDM(:,(stop+1):end) = [];


end


