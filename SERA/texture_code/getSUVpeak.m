function [SUVpeak] = getSUVpeak(RawImg,ROI,pixelW,sliceTh)
% -------------------------------------------------------------------------
% function [SUVpeak] = getSUVpeak(RawImg,ROI,pixelW,sliceS)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the SUVpeak or local intensity peak in two ways:
% Local -> sphere around the maximum voxel
% Global -> sphere around each voxel
%
% --> This function is compatible with 2D analysis
% -------------------------------------------------------------------------
% INPUTS:
% - RawImg: This is the original ROI. This code will further trim it to the
%        ROI region. But it preserves the non-ROI values required for
%        SUVpeak calculation.
% - ROI: The whole ROI. It will be trimmed inside this code.
% - pixelW and SliceS: pixel width and slice thickness
%
% -------------------------------------------------------------------------
% OUTPUTS:
% - SUVpeak, a 2 element array consisting local and global SUVpeak.
% -------------------------------------------------------------------------
% AUTHOR(S):
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2017
% - Revision: March 2019
% -------------------------------------------------------------------------
% --> Copyright (c) 2007-2017, Saeed Ashrafinia, Johns Hopkins University
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

% [boxBound] = computeBoundingBox(ROI);
% ROIBox = ROI(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
% ImgBox = RawImg(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
% ImgRawROI = RawImg .* ROI;

%% Create Spheroid
% SUVpeak (based on the definition of local intensity peak)
% based on definition of 1cm3 sephere
R = ((3/(4*pi))^(1/3) * 10 ./ [pixelW,sliceTh]);
SPH=zeros(2*floor(R(1))+1,2*floor(R(1))+1,2*floor(R(2))+1);

%% Create the sphere
% This newer method insures the distance from the "center of voxels" are
% within the required threshold from the max voxel.
[x,y,z] = meshgrid(pixelW.*((-ceil(size(SPH,2)/2)+0.5:floor(size(SPH,2)/2))-0.5),...
                   pixelW.*((-ceil(size(SPH,1)/2)+0.5:floor(size(SPH,1)/2))-0.5),...
                   sliceTh.*((-ceil(size(SPH,3)/2)+0.5:floor(size(SPH,3)/2))-0.5)); % In mm
tmpsph = sqrt((x-x(1,ceil(size(x,1)/2))).^2 + (y-y(ceil(size(y,2)/2),1)).^2 + (z-z(1,1,ceil(size(z,3)/2))).^2);
tmpsph(find(tmpsph > ((3/(4*pi))^(1/3)*10))) = NaN; %#ok<*FNDSB>
SPH = tmpsph;
SPH(find(~isnan(tmpsph))) = 1; 

%% New method based on Alex's suggestion, much faster
R = floor(R);
ImgRawROIpadded = padarray(RawImg,[R(1) R(1) R(2)],NaN); % Pad the image with extra NaNs to ensure we cover spheres centered at voxels in the edge of the image
ImgRawROIpadded(find(isnan(ImgRawROIpadded)))=0;
SPH(find(isnan(SPH)))=0;
% Perform convolution instead of for loop
C=convn(ImgRawROIpadded, SPH./sum(SPH(:),'omitnan'),'valid');

T=[RawImg(:), ROI(:), C(:)];
T1=T;

T1(find(isnan(T1(:,1))),:)=[];
T2=T1(find(T1(:,2)~=0),:);
[~,maxind]=max(T2(:,1));
SUVpeak = max(T2(maxind,3));    % Local intensity peak
SUVpeak(2,1) = max(T2(:,3));      % Global intensity peak




% % This is the old method
% %% Local Intensity Peak
% ROInan = ROI; % Create an ROI
% ROInan(ROI==0) = NaN;
% R = floor(R);
% ImgRawROIpadded = padarray(RawImg,[R(1) R(1) R(2)],NaN); % Pad the image with extra NaNs to ensure we cover spheres centered at voxels in the edge of the image
% ImgROIpadded = ImgRawROIpadded .* padarray(cast(ROInan,class(ImgRawROIpadded)),[R(1) R(1) R(2)],NaN); % Pad the ROI too, keep the class type
% MaxIntensity = max(ImgROIpadded(:)); % Find the maximum value of voxels
% MaxIndx = find(ImgROIpadded == MaxIntensity); % Find voxels with this max value
% [indMaxX,indMaxY,indMaxZ] = ind2sub(size(ImgROIpadded),MaxIndx); % Get their x,y,z
% peakTmpLocal = []; % to store peak values corresponding to all voxels that have max value
% for ind = 1:length(MaxIndx)
%     peakVol = ImgRawROIpadded(indMaxX(ind)-R(1):indMaxX(ind)+R(1) , indMaxY(ind)-R(1):indMaxY(ind)+R(1) , indMaxZ(ind)-R(2):indMaxZ(ind)+R(2)) .* SPH;
%     peakTmpLocal = cat(2,peakTmpLocal , mean(peakVol(:),'omitnan'));
% end
% p1 = max(peakTmpLocal); % Save the maximum peak value as the local intensity peak.
% 
% 
% %% Global Intensity Peak
% if isclcGlobPeak
%     MaxIndx = find(~isnan(ImgROIpadded));
%     [indMaxX,indMaxY,indMaxZ] = ind2sub(size(ImgROIpadded),MaxIndx);
%     peakTmpGlob = [];
%     peakTmpLocal = [];
%     for ind = 1:length(MaxIndx)
%         peakVol = ImgRawROIpadded(indMaxX(ind)-R(1):indMaxX(ind)+R(1) , indMaxY(ind)-R(1):indMaxY(ind)+R(1) , indMaxZ(ind)-R(2):indMaxZ(ind)+R(2)) .* SPH;
%         peakTmpGlob = cat(2,peakTmpGlob , mean(peakVol(:),'omitnan'));
%         if ImgRawROIpadded(indMaxX(ind) , indMaxY(ind) , indMaxZ(ind)) == MaxIntensity
%             peakTmpLocal = cat(2,peakTmpLocal , mean(peakVol(:),'omitnan'));
%         end
%     end
%     SUVpeak = max(peakTmpLocal);
% %     if p1 ~=SUVpeak, error('SUVpeak local has an issue'); end
% else
%     peakTmpGlob = p1; 
% end
% SUVpeak = cat(1, p1 , max(peakTmpGlob));
% 

















