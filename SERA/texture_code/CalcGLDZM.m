function [ArrayOut] = CalcGLDZM(GLDZM,ROI)
% -------------------------------------------------------------------------
% function [SM_f, SS_f] = getGLCM2Dtex(ROIonly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates 16 GLDZM features for both 2D and 3D. 
% -------------------------------------------------------------------------
% INPUTS:
% - GLDZM: Grey Level Distant Zone Matrix
% -------------------------------------------------------------------------
% OUTPUTS:
% - FeatMatrixout: an array of 16 GLDZM fetures
% -------------------------------------------------------------------------
% AUTHOR(S):    
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
% - Revision: July 2017
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

nG = size(GLDZM,1);             % Number of Gray Levels
nD = size(GLDZM,2);             % Number of dependence + 1 (j = k+1)
nS = sum(GLDZM(:));             % The number of neighbourhoods.
Di = sum(GLDZM , 2);            % The number of neighbourhoods with discretised grey level i, essentially constituting a grey level histogram of the volume or slice
Dj = sum(GLDZM , 1);            % The number of neighbourhoods with dependence j, regardless of grey level.
nV = numel(find(~isnan(ROI)));  % Number of all voxels in the ROI
% FeatMatrixout = [];
% for k = 1:size(GLDZM,3)
%     tmpGLDZM = squeeze(GLDZM(:,:,k));
%     nS = sum(tmpGLDZM(:));             % Sum of all elements in GLDZM
% %     nV = numel(find(~isnan(ROI(:,:,k))));  % Number of all voxels in the ROI
%     FeatMatrixout = cat(2,FeatMatrixout , GLDZMfeatHandler(tmpGLDZM,nG,nD,nS,nV)');
% end  
% end
% 
% 
% function [ArrayOut] = GLDZMfeatHandler(GLDZM,nG,nD,nS,nV)

% Begin calling functions
ArrayOut(1)     = SmallDistEmph(Dj,nD,nS);  %N
ArrayOut(2)     = LargeDistEmph(Dj,nD,nS);  %N
ArrayOut(3)     = LowGLCountRmph(Di,nG,nS);  %*
ArrayOut(4)     = HighGLCountEmph(Di,nG,nS);  %*
ArrayOut(5)     = SmallDistLowGLEmph(GLDZM,nG,nD,nS);  %N
ArrayOut(6)     = SmallDistHighRL(GLDZM,nG,nD,nS);  %N
ArrayOut(7)     = LargeDistLowGLEmph(GLDZM,nG,nD,nS);  %N
ArrayOut(8)     = LargeDistHighGLEmph(GLDZM,nG,nD,nS);  %N
ArrayOut(9)     = GLnonUnif(Di,nS);  %*
ArrayOut(10)    = GLnonUnifNormzd(Di,nS);  %*
ArrayOut(11)    = ZoneDistNonUnif(Dj,nS);  %N
ArrayOut(12)    = ZoneDistNonUnifNormzd(Dj,nS);  %N
ArrayOut(13)    = ZonePercentage(nS,nV);  %*
ArrayOut(14)    = GLVar(GLDZM,nG,nD,nS);  %*
ArrayOut(15)    = ZoneDistVar(GLDZM,nG,nD,nS);  %N
ArrayOut(16)    = ZoneDistEntropy(GLDZM,nS);  %




end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature Calculation functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (SDE) Small distance emphasis (1)
function [f_SDE] = SmallDistEmph(Dj,nD,nS)
% This feature emphasises short run lengths

J = 1:nD;
tmp = Dj ./ J.^2;
f_SDE = sum(tmp)/nS;
end

%% (LDE) Large distance emphasis (2)
function [f_LDE] = LargeDistEmph(Dj,nD,nS)
% This feature emphasises long run lengths

J = 1:nD;
tmp = Dj .* J.^2;
f_LDE = sum(tmp)/nS;
end

%% (LGLZE) Low grey level zone emphasis (3)
function [f_LGLZE] = LowGLCountRmph(Di,nG,nS)
% This feature is a grey level analogue to f_SRE (Chu et al., 1990).
% Instead of low run lengths, low grey levels are emphasised

I = (1:nG)';
tmp = Di ./ I.^2;
f_LGLZE = sum(tmp)/nS;
end

%% (HGLZE) High grey level zone emphasis (4)
function [f_HGLZE] = HighGLCountEmph(Di,nG,nS)
% The high grey level run emphasis feature is a grey level analogue to
% Frlm:lre (Chu et al., 1990). The feature emphasises high grey levels.

I = (1:nG)';
tmp = Di .* I.^2;
f_HGLZE = sum(tmp)/nS;
end

%% (SDLGLE) Small distance low grey level emphasis (5)
function [f_SRLGLE] = SmallDistLowGLEmph(GLDZM,nG,nD,nS)
% This feature emphasises runs in the upper left quadrant of the GLDZM,
% where short run lengths and low grey levels are located

[I , J] = ndgrid(1:nG , 1:nD);
tmp = GLDZM ./ (I.^2 .* J.^2);
f_SRLGLE = sum(tmp(:))/nS;
end

%% (SDHGLE) Small distance high grey level emphasis (6)
function [f_SRHGLE] = SmallDistHighRL(GLDZM,nG,nD,nS)
% This feature emphasises runs in the lower left quadrant of the GLDZM,
% where short run lengths and high grey levels are located

[I , J] = ndgrid(1:nG , 1:nD);
tmp = I.^2 .* GLDZM ./ J.^2;
f_SRHGLE = sum(tmp(:))/nS;
end

%% (LDLGLE) Large distance low grey level emphasis (7)
function [f_LRLGLE] = LargeDistLowGLEmph(GLDZM,nG,nD,nS)
% This feature emphasises runs in the upper right quadrant of the GLDZM,
% where long run lengths and low grey levels are located 

[I , J] = ndgrid(1:nG , 1:nD);
tmp = J.^2 .* GLDZM ./ I.^2;
f_LRLGLE = sum(tmp(:))/nS;
end

%% (LDHGLE) Large distance high grey level emphasis (8)
function [f_LRHGLE] = LargeDistHighGLEmph(GLDZM,nG,nD,nS)
% This feature emphasises runs in the lower right quadrant of the GLDZM,
% where long run lengths and high grey levels are located  

[I , J] = ndgrid(1:nG , 1:nD);
tmp = (J.^2) .* (I.^2) .* GLDZM;
f_LRHGLE = sum(tmp(:))/nS;
end

%% (GLNU) Grey level non-uniformity (9)
function [f_GLNU] = GLnonUnif(Di,nS)
% This feature assesses the distribution of runs over the grey values
% (Galloway, 1975). The feature value is low when runs are equally
% distributed along grey levels. 

f_GLNU = sum(Di.^2) / nS;
end

%% (GLNUN) Grey level non-uniformity normalized (10)
function [f_GLNUN] = GLnonUnifNormzd(Di,nS)
% This is a normalised version of the grey level non-uniformity feature.

f_GLNUN = sum(Di.^2) / nS.^2;
end

%% (ZDNU) Zone distance non-uniformity (11)
function [f_ZDNU] = ZoneDistNonUnif(Dj,nS)
% This features assesses the distribution of runs over the run lengths
% (Galloway, 1975). The feature value is low when runs are equally
% distributed along run lengths. 

f_ZDNU = sum(Dj.^2) / nS;
end

%% (ZDNUN) Zone distance non-uniformity normalised (12)
function [f_ZDNUN] = ZoneDistNonUnifNormzd(Dj,nS)
% This is normalised version of the run length non-uniformity feature.

f_ZDNUN = sum(Dj.^2) / nS.^2;
end

%% (ZP) Zone percentage (13)
function [f_ZP] = ZonePercentage(nS,nV)
% This feature assesses the fraction of the number of realised runs and the
% maximum number of potential runs (Galloway, 1975). Strongly linear or
% highly uniform ROI volumes produce a low run percentage.  

f_ZP = nS / nV;
end

%% (GLV) Grey level variance (14)
function [f_GLV] = GLVar(GLDZM,nG,nD,nS)
% This feature estimates the variance in runs for the grey levels.

Pij = GLDZM / nS;
[I , ~] = ndgrid(1:nG , 1:nD);
mu = sum(sum(I .* Pij));
tmp = (I - mu).^2 .* Pij;
f_GLV = sum(tmp(:));
end

%% (ZDV) Zone distance variance (15)
function [f_ZDV] = ZoneDistVar(GLDZM,nG,nD,nS)
% This feature estimates the variance in runs for run lengths.

Pij = GLDZM / nS;
[~ , J] = ndgrid(1:nG , 1:nD);
mu = sum(sum(J .* Pij));
tmp = (J - mu).^2 .* Pij;
f_ZDV = sum(tmp(:));
end

%% (ZDE) Zone distance entropy (16)
function [f_ZD] = ZoneDistEntropy(GLDZM,nS)
% It calculates the run entropy.

Pij = GLDZM / nS;
tmp = Pij .* log2(Pij + realmin);
f_ZD = -sum(tmp(:));
end



