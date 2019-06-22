function [VDs, K] = getVDs(ROIbox,ROIboxV,pixelW,sliceS,MV,Surface )
% -------------------------------------------------------------------------
% function [VDs, K] = getVDs(ROIbox,ROIboxV,pixelW,sliceS,MV,Surface )
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function several morpholical features for an ROI, including major,
% minor and least axis lengths; flatness and elongation, various volume and
% area densities as defined in ISBI documentation. It also calculates the
% vertices of convex Hull of the ROI.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIbox: The smallest box containing the 3D morphological ROI. Voxels
%            outside the ROI are set to NaNs.
% - ROIboxV: A Nx3 matrix of vertices of ROI. Use MarchingCubes function
%            to generate ROIboxV.
% - pixelW: width of the voxel in the X (=Y) direction
% - sliceS: Slice thickness of the voxel in the Z direction
% - MV: tumor volume
% - Surface: tumor surface
% -------------------------------------------------------------------------
% OUTPUTS:
% - VDs: An array of calculated volume densities as formulated below. 
% - K: vertices of the Convex Hull of the ROI
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

%% Find eigenvalues
[x,y,z]=ind2sub(size(ROIbox),find(~isnan(ROIbox)));
% [x,y,z]=ind2sub(size(ROIbox),find(ROIbox~=0));
% [~,~,EIGs]=pca([x*pixelW, y*pixelW, z*sliceS],'Centered','off');
EIGs=eig(cov([x*pixelW, y*pixelW, z*sliceS]));
% EIGs = diag(getPCAeig(ROIbox, pixelW,sliceS))'; 
EIGs = flip(EIGs)';
if numel(EIGs) == 1
    EIGs = repmat(EIGs,1,3);
end
L1 = 2*sqrt(EIGs(1));  L2 = 2*sqrt(EIGs(2));  L3 = 2*sqrt(EIGs(3)); 

%% calculate major, minor, least axis, elongation, flatness
MALs = 4*sqrt(EIGs);                    % These are Major, Milor and Least axis length
Elong = sqrt(EIGs(2) / EIGs(1));        % Elongation
Fltns = sqrt(EIGs(3) / EIGs(1));        % Flatness

%% calculate volume and area densities
[nx, ny, nz] = size(ROIbox);
% Volume density - axis-aligned bounding box
VDaabb = MV / prod([nx ny nz] .* [pixelW,pixelW,sliceS]);      
% Area density - axis-aligned bounding box
ADaabb = Surface / (2*nx*pixelW*ny*pixelW + 2*ny*pixelW*nz*sliceS + 2*nx*pixelW*nz*sliceS);  
% Volume density - approximate enclosing ellipsoid
VDaee  = 3*MV / (4*pi*2^3*prod(sqrt(EIGs)));                    

%% calculate area density - approximate enclosing ellipsoid:
alpha = sqrt(1-L2^2/L1^2);
beta = sqrt(1-L3^2/L1^2);
MaxV = 20;
v0 = 0:MaxV;
sumtmp = []; 
for k = 0: MaxV
    sumtmp = cat(2,sumtmp , (alpha*beta).^k ./ (1-4*k.^2) .* legendre(k , min(1,(alpha^2+beta^2)/(2*alpha*beta)))');
end
aee = 4*pi*L1*L2 * sum(sumtmp);
ADaee = Surface / aee;

%% Volume density - minimum volume enclosing ellipsoid
% [~,EIGs]=pcafun(ROIboxV); 
[~,EIGs]=eig(cov(ROIboxV)); 
EIGs = sum(EIGs , 1); 
L1 = 2*sqrt(EIGs(1));  L2 = 2*sqrt(EIGs(2));  L3 = 2*sqrt(EIGs(3)); 
VDmvee = MV / (4*pi*L1*L2*L3/3);

%% calculate area density - minimum volume enclosing ellipsoid:
alpha = sqrt(1-L2^2/L1^2);
beta = sqrt(1-L3^2/L1^2);
try
    aee = 4*pi*L1*L2 * sum((alpha*beta).^v0 ./ (1-4*v0.^2) .* legendre(MaxV , min(1,(alpha^2+beta^2)/(2*alpha*beta)))');
catch
    EIGs = flip(EIGs)';
    
    L1 = 2*sqrt(EIGs(1));  L2 = 2*sqrt(EIGs(2));  L3 = 2*sqrt(EIGs(3));
    VDmvee = MV / (4*pi*L1*L2*L3/3);
    alpha = sqrt(1-L2^2/L1^2);
    beta = sqrt(1-L3^2/L1^2);
    aee = 4*pi*L1*L2 * sum((alpha*beta).^v0 ./ (1-4*v0.^2) .* legendre(MaxV , min(1,(alpha^2+beta^2)/(2*alpha*beta)))');
end
    
ADmvee = Surface / aee;

%% Oriented minimum bounding box
try
[~,~,Vombb,Aombb,~] = minboundbox(x*pixelW, y*pixelW, z*sliceS,'volume',1);
catch
    disp('This is a 2D or 1D ROI. Switch to 2D Convex Hull and Bounding Box calculation.')
    [xtmp,ytmp]=ind2sub(size(squeeze(ROIbox)),find(~isnan(squeeze(ROIbox))));
    try
        if size(ROIbox,3) == 1
            [~ , Aombb, Vombb] = minBoundingBox2D([xtmp , ytmp]',pixelW,pixelW,sliceS);
        elseif size(ROIbox,2) == 1
            [~ , Aombb, Vombb] = minBoundingBox2D([xtmp , ytmp]',pixelW,sliceS,pixelW);
        elseif size(ROIbox,1) == 1
            [~ , Aombb, Vombb] = minBoundingBox2D([xtmp , ytmp]',sliceS,pixelW,pixelW);
        else
            warning('Min bounding box does not respond (no Convex Hull available). Set OMBB = AABB....')
            Vombb = prod([nx ny nz] .* [pixelW,pixelW,sliceS]);
            Aombb = (2*nx*pixelW*ny*pixelW + 2*ny*pixelW*nz*sliceS + 2*nx*pixelW*nz*sliceS);
            
        end
    catch
        disp('Min bounding box does not respond (probably a 1D ROI). Set OMBB = AABB....')
        Vombb = prod([nx ny nz] .* [pixelW,pixelW,sliceS]);
        Aombb = (2*nx*pixelW*ny*pixelW + 2*ny*pixelW*nz*sliceS + 2*nx*pixelW*nz*sliceS);
    end
        
end
VDombb = MV / Vombb;
ADombb = Surface / Aombb;

%% calculate area and volume density of Convex Hull:
[K,Vch] = convhulln(ROIboxV);%[x.*pixelW y.*pixelW z.*sliceS],{'Qt'}); 
Solidity = MV / Vch; 
CHarea = sum(sqrt(sum(( ...
[ROIboxV(K(:,1),2).*ROIboxV(K(:,2),3) - ROIboxV(K(:,1),3).*ROIboxV(K(:,2),2) ...
ROIboxV(K(:,1),3).*ROIboxV(K(:,2),1) - ROIboxV(K(:,1),1).*ROIboxV(K(:,2),3)  ...
ROIboxV(K(:,1),1).*ROIboxV(K(:,2),2) - ROIboxV(K(:,1),2).*ROIboxV(K(:,2),1)] + ...
[ROIboxV(K(:,2),2).*ROIboxV(K(:,3),3) - ROIboxV(K(:,2),3).*ROIboxV(K(:,3),2) ...
ROIboxV(K(:,2),3).*ROIboxV(K(:,3),1) - ROIboxV(K(:,2),1).*ROIboxV(K(:,3),3)  ...
ROIboxV(K(:,2),1).*ROIboxV(K(:,3),2) - ROIboxV(K(:,2),2).*ROIboxV(K(:,3),1)] + ...
[ROIboxV(K(:,3),2).*ROIboxV(K(:,1),3) - ROIboxV(K(:,3),3).*ROIboxV(K(:,1),2) ...
ROIboxV(K(:,3),3).*ROIboxV(K(:,1),1) - ROIboxV(K(:,3),1).*ROIboxV(K(:,1),3)  ...
ROIboxV(K(:,3),1).*ROIboxV(K(:,1),2) - ROIboxV(K(:,3),2).*ROIboxV(K(:,1),1)]).^2,2))) /2;
CH_AD = Surface/CHarea;


%% Export all 
VDs = [MALs, Elong ,Fltns ,VDaabb ,ADaabb ,VDombb, ADombb, VDaee ,ADaee, VDmvee, ADmvee, Solidity, CH_AD]';
