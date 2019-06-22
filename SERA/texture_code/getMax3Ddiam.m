function Max3Ddiam = getMax3Ddiam(ROIboxV, ConvHullV)
% -------------------------------------------------------------------------
% function Max3Ddiam = getMax3Ddiam(ROIboxV, ConvHullV)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes maximum 3d diameter of an ROI. For ROIs with small
% number of voxels (<5000) it uses the exact ROI, and for larger ROIs it
% uses the ROI from convex Hull. 
% -------------------------------------------------------------------------
% INPUTS:
% - ROIboxV: A Nx3 matrix of vertices of ROI. Use MarchingCubes function
%            to generate ROIboxV.
% - ConvHullV: A Nx3 matrix of vertices of the convex hull of the ROIs. Use
%              convhulln function to generate ConvHullV.
% -------------------------------------------------------------------------
% OUTPUTS:
% - Max3Ddiam: maximum 3d diameter of an ROI 
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

% This function calculates the maximum diameter in 3D of an ROI.
% Written in two methods. The fastest is chosen.
try
NumVer          = size(ROIboxV,1);
ROIboxV_single  = single(ROIboxV);
r2              = ROIboxV_single';
r2              = reshape((repmat(ROIboxV_single,1,NumVer) - repmat(r2(:)',NumVer,1)).^2 , [NumVer,3,NumVer]);
rss             = sqrt(squeeze(sum(r2,2)));
Max3Ddiam       = max(rss(:));

catch
    warning('Problem with calculating Max 3D diameter, most probably issue with "low memory". An alternative "approximate" way using Convex Hull will be selected.');
    
    NumVer = size(ConvHullV,1);
    r1=[ROIboxV(ConvHullV(:,1),1),ROIboxV(ConvHullV(:,2),2),ROIboxV(ConvHullV(:,3),3)];
    r1m = repmat(r1,1,NumVer);
    r2=r1';
%     r2=r2(:)';
    r2m=repmat(r2(:)',NumVer,1);
    rm=(r1m-r2m).^2;
    rs=reshape(rm,[],3);
    rss=sqrt(sum(rs,2));
    Max3Ddiam = max(rss);
end

% NumVer = size(ROIboxV,1);
% if NumVer < 5000
%     % Method 1
%     r1m = repmat(ROIboxV,1,NumVer);
%     r2=ROIboxV';
% %     r2=r2(:)';
%     r2=repmat(r2(:)',size(ROIboxV,1),1);
%     rm=(r1m-r2).^2;
%     rs=reshape(rm,[NumVer,3,NumVer]);
%     rss=sqrt(squeeze(sum(rs,2)));
%     Max3Ddiam = max(rss(:));
% else
%     NumVer = size(ConvHullV,1);
%     r1=[ROIboxV(ConvHullV(:,1),1),ROIboxV(ConvHullV(:,2),2),ROIboxV(ConvHullV(:,3),3)];
%     r1m = repmat(r1,1,NumVer);
%     r2=r1';
% %     r2=r2(:)';
%     r2m=repmat(r2(:)',NumVer,1);
%     rm=(r1m-r2m).^2;
%     rs=reshape(rm,[],3);
%     rss=sqrt(sum(rs,2));
%     Max3Ddiam = max(rss);
% end

% % % Method 2
% % c=mat2cell(ROIboxV,[linspace(1,1,size(ROIboxV,1))],[3]);
% % rv = repmat(c,1,size(ROIboxV,1));
% % rvt = repmat(c',size(ROIboxV,1),1);
% % s=cellfun(@norm,cellfun(@minus,rv,rvt,'UniformOutput',false),'UniformOutput',false);
% % max(cell2mat(s(:)))

end