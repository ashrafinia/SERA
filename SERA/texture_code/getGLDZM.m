function [GLDZM] = getGLDZM(ROIBox,ROIMask,levels)
% -------------------------------------------------------------------------
% function [GLDZM] = getGLDZM(ROIonly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes the Neighbouring grey level dependence Matrix
% (GLDZM) of the region of interest (ROI) of an input volume. The input
% volume is assumed to be isotropically resampled. Only one GLDZM is
% computed per scan.    
%
% --> This function is compatible with 2D analysis 
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
% - GLDZM: Gray-Level Run-Length Matrix of 'ROIonly'.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Saeed Ashrafinia 
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2017
% - Revision: Mar 2019
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

% PRELIMINARY
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
levelTemp = max(levels) + 1;
ROIBox(isnan(ROIBox)) = levelTemp;
levels = [levels,levelTemp];


% QUANTIZATION EFFECTS CORRECTION
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
uniqueVect = round(levels*adjust)/adjust;
ROIBox = round(ROIBox*adjust)/adjust;
% ROIOnlyP = padarray(ROIOnly,[1 1 1],'both');
nGL = length(levels) - 1;


% INITIALIZATION
% nInit = numel(ROIBox);
GLDZM = zeros(nGL,ceil(max(size(ROIBox)/2)));

unqGLs = unique(ROIBox);       % find unique GLs within the current ROI
unqGLs(unqGLs==levelTemp) = []; % exclude the non-ROI voxels from unqGLs

% COMPUTATION OF GLDZM
if numel(size(ROIMask)) == 3
    ROIOnlyP = padarray(ROIMask,[1 1 1],'both');
    ROIBoxP = padarray(ROIBox,[1 1 1],'both');
    ROIBoxP(ROIBoxP == 0) = levelTemp;
    
    % Create the distance map
    % find distance of all voxels to any edge or any non-ROI voxel 
    try % if ROI type is not single or double, apparetnly returns error. so is we get error, we set the class to single.
    distmap = bwdistNew(ROIOnlyP);
    catch
    distmap = bwdistNew(single(ROIOnlyP));
    end
    distmap(distmap==0)=NaN;
    
    for gl = 1:numel(unqGLs)
        % Create the zone size matrix
        temp = ROIOnlyP;
        temp(ROIBoxP~=unqGLs(gl)) = 0;
        temp(ROIBoxP==unqGLs(gl)) = 1;
        connObjects = bwconncomp(temp,26);  % find the connected components
        nConn = length(connObjects.PixelIdxList);
        
        for c = 1:nConn
            
            tmpROI = zeros(size(distmap));
            tmpROI(connObjects.PixelIdxList{c}) = 1;
            tmpROI(tmpROI==0) = NaN;
            Dist = min(min(min(tmpROI .* distmap,[],'omitnan')));
            idx = find(unqGLs(gl)==uniqueVect,1,'first');
            try
                GLDZM(idx,Dist) = GLDZM(idx,Dist) + 1;
            catch
                disp(['idx = ',num2str(idx)]);
                disp(['gl = ',num2str(gl)]);
                disp(['idx = ',num2str(uniqueVect(idx))]);
                disp(['unqGLs(gl) = ',num2str(unqGLs(gl))]);
            end
        end
    end
    
elseif numel(size(ROIMask)) <= 2
    ROIOnlyP = padarray(ROIMask,[1 1],'both');
    ROIBoxP = padarray(ROIBox,[1 1],'both');
    ROIBoxP(ROIBoxP == 0) = levelTemp;
    
    % Create the distance map
    try % if ROI type is not single or double, apparetnly returns error. so is we get error, we set the class to single.
    distmap = bwdistNew(ROIOnlyP);
    catch
    distmap = bwdistNew(single(ROIOnlyP));
    end
    distmap(distmap==0)=NaN;
    
    for gl = 1:numel(unqGLs)
        % Create the zone size matrix
        temp = ROIOnlyP;
        temp(ROIBoxP~=unqGLs(gl)) = 0;
        temp(ROIBoxP==unqGLs(gl)) = 1;
        connObjects = bwconncomp(temp,8);
        nConn = length(connObjects.PixelIdxList);
        
        for c = 1:nConn
            tmpROI = zeros(size(distmap));
            tmpROI(connObjects.PixelIdxList{c}) = 1;
            tmpROI(tmpROI==0) = NaN;
            Dist = min(min(tmpROI .* distmap,[],'omitnan'));
            idx = find(unqGLs(gl)==uniqueVect,1,'first');
            
            try
                GLDZM(idx,Dist) = GLDZM(idx,Dist) + 1;
            catch
                disp(['idx = ',num2str(idx)]);
                disp(['gl = ',num2str(gl)]);
                disp(['idx = ',num2str(uniqueVect(idx))]);
                disp(['unqGLs(gl) = ',num2str(unqGLs(gl))]);
            end
        end
    end
end


% REMOVE UNECESSARY COLUMNS
stop = find(sum(GLDZM),1,'last');
GLDZM(:,(stop+2):end) = [];






end



% % COMPUTATION OF GLDZM
% if numel(size(ROIOnly))==3
%     temp = ROIOnly;
%     [I,J,K] = ind2sub(size(ROIOnly),find(~isnan(ROIOnly)));
%     for gl = 1:nGL
%         temp(ROIOnly~=uniqueVect(gl)) = 0;
%         temp(ROIOnly==uniqueVect(gl)) = 1;
%         connObjects = bwconncomp(temp,26);
%         %     numPix = cellfun(@numel,connObjects.PixelIdxList);
%         %     GLDZM(i,numPix) = GLDZM(i,numPix) + 1;
%         
%         nConn = length(connObjects.PixelIdxList);
%         for c = 1:nConn
%             [zi,zj,zk] = ind2sub(size(ROIOnly),cell2mat(connObjects.PixelIdxList(c)));
%             Dist = min(min([zi-1, I(end)-zi, zj-1, J(end)-zj, zk-1, K(end)-zk]));        % Find the minimum distance to one of the edges
%             GLDZM(gl,Dist+1) = GLDZM(gl,Dist+1) + 1;
%         end
%     end
%     
% elseif numel(size(ROIOnly)) <= 2
%     temp = ROIOnly;
%     [I,J] = ind2sub(size(ROIOnly),find(~isnan(ROIOnly)));
%     for gl = 1:nGL
%         temp(ROIOnly~=uniqueVect(gl)) = 0;
%         temp(ROIOnly==uniqueVect(gl)) = 1;
%         connObjects = bwconncomp(temp,26);
%         nConn = length(connObjects.PixelIdxList);
%         for c = 1:nConn
%             [zi,zj] = ind2sub(size(ROIOnly),cell2mat(connObjects.PixelIdxList(c)));
%             Dist = min(min([zi-1, I(end)-zi, zj-1, J(end)-zj]));        % Find the minimum distance to one of the edges
%             GLDZM(gl,Dist+1) = GLDZM(gl,Dist+1) + 1;
%         end
%     end
% end
% 
% 
