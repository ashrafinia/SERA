function [MoranI, GearyC]= getAutoCorrs(ROIints , xI,yI,zI , pixelW,sliceS , Mean)
% -------------------------------------------------------------------------
%% % function MoranI = getMoranI(ROIints , Mean)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates the Moran's I index and Geary's C measure for
% the intensity ROI. 
% Moran's I index is an indicator of spatial autocorrelation (Moran, 1950).
% Geary's C measures spatial autocorrelation, like Moran's I index (Geary,
% 1954). Geary's C however, directly measures grey level differences
% between voxels and is more sensitive to local spatial autocorrelation.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIints: The smallest box containing the 3D intensity ROI with
%             their intensities. Voxels outside the ROI are set to NaNs.
% - Mean: Mean of the intensity box
% -------------------------------------------------------------------------
% OUTPUTS:
% - MoranI: Moran's I index
% - GearyC: Geary's C measure 
% -------------------------------------------------------------------------
% AUTHOR(S): Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: Nov 2017
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

Xgl = single(ROIints(~isnan(ROIints)));
Points = single([xI*pixelW,yI*pixelW,zI*sliceS]);
%taking care of 1D
if size(Xgl,1) == 1
    try
    Xgl = Xgl';
    catch
        Xgl = squeeze(Xgl);
    end
end
if size(Points,1) == 1
    Points = single([xI'*pixelW,yI'*pixelW,zI'*sliceS]);
end
nV = length(Xgl);
MNum = 0; GNum = 0; sumW = 0;

for r = 1:nV
    W = 1./ sqrt((Points(:,1)-Points(r,1)).^2 + (Points(:,2)-Points(r,2)).^2 + (Points(:,3)-Points(r,3)).^2);
    W(r) = 0;
    try % taking care of 1D/2D errors
        MNum = MNum + W' * ((Xgl(r)-Mean) .* (Xgl-Mean));
    catch
        MNum = MNum + squeeze(W)' * squeeze((Xgl(r)-Mean) .* (Xgl-Mean));
    end
    GNum = GNum + sum(squeeze(W) .* squeeze((Xgl(r) - Xgl) .^ 2));
    sumW = sumW + sum(W);
end
denom = sum((Xgl - Mean).^2);
MoranI = nV/sumW * MNum / denom; 
GearyC = (nV-1)/(2*sumW) * GNum / denom;
end

% %%
% Xgl = single(ROIints(~isnan(ROIints)));
% Points = single([xI*pixelW,yI*pixelW,zI*sliceS]);
% nV = length(Xgl);
% nV_4 = floor(nV/4); 
% MNum = 0; GNum = 0; sumW = 0; 
% 
% %  = [{1:nV_4}, {(nV_4+1):(2*nV_4)}, {(2*nV_4+1):(3*nV_4)}, {(3*nV_4+1):nV}];
% for mx = 1:4
%     if mx~=4, IndX = ((mx-1)*nV_4+1) : (mx*nV_4); else  IndX = ((mx-1)*nV_4+1) : nV; end
%     for my = 1:4
%         if my~=4, IndY = ((my-1)*nV_4+1) : (my*nV_4); else  IndY = ((my-1)*nV_4+1) :nV; end 
%         try
%         W = sqrt((repmat(Points(IndX,1),1,numel(IndY)) - repmat(Points(IndY,1)',numel(IndX),1)).^2 + ...
%                  (repmat(Points(IndX,2),1,numel(IndY)) - repmat(Points(IndY,2)',numel(IndX),1)).^2 + ...
%                  (repmat(Points(IndX,3),1,numel(IndY)) - repmat(Points(IndY,3)',numel(IndX),1)).^2);
%         catch
%             disp('Wrong'); 
%         end
%         W = 1 ./ W;
%         if mx == my
%             W(1== eye(size(W))) = 0;
%         end
%         MNum = MNum + sum(sum(W .* repmat(Xgl(IndX) - Mean , 1,numel(IndY)) .* repmat((Xgl(IndY)-Mean)' , numel(IndX),1)));
%         GNum = GNum + sum(sum(W .* (repmat(Xgl(IndX), 1,numel(IndY)) - repmat(Xgl(IndY)' , numel(IndX),1)).^2));
%         sumW = sumW + sum(W(:));
%     end
% end
% 
% denom = sum((Xgl - Mean).^2);
% MoranI = nV/sumW * MNum / denom; 
% GearyC = (nV-1)/(2*sumW) * GNum / denom;
% 
% 
% %%
% Xgl = single(ROIints(~isnan(ROIints)));
% Points = single([xI*pixelW,yI*pixelW,zI*sliceS]);
% nV = length(Xgl);
% 
% W = 1 ./ sqrt((Points(2:end,1)-Points(1,1)).^2 + (Points(2:end,2)-Points(1,2)).^2 + (Points(2:end,3)-Points(1,3)).^2);
% MNum = sum((Xgl(1)-Mean) * W .* (Xgl(2:end)-Mean));
% GNum = sum(W .* (Xgl(1) - Xgl(2:end)).^2);
% sumW = sum(W);
%     
% for j = 2:nV-1
%     interval = [1:j-1 , j+1:nV];
%     W = 1 ./ sqrt((Points(interval,1)-Points(j,1)).^2 + (Points(interval,2)-Points(j,2)).^2 + (Points(interval,3)-Points(j,3)).^2 );
%     MNum = MNum + sum((Xgl(j)-Mean) * W .* (Xgl(interval)-Mean));
%     GNum = GNum + sum(W .* (Xgl(j) - Xgl(interval)).^2);
%     sumW = sumW + sum(W);
%     
% end
% 
% W = 1 ./ sqrt((Points(1:end-1,1)-Points(end,1)).^2 + (Points(1:end-1,2)-Points(end,2)).^2 + (Points(1:end-1,3)-Points(end,3)).^2);
% MNum = MNum + sum((Xgl(end)-Mean) * sumW .* (Xgl(1:end-1)-Mean));
% GNum = GNum + sum(sumW .* (Xgl(end) - Xgl(1:end-1)).^2);
% sumW = sumW + sum(W);
% 
% 
% denom = sum((Xgl - Mean).^2);
% MoranI = nV/sumW * MNum / denom; 
% GearyC = (nV-1)/(2*sumW) * GNum / denom;
