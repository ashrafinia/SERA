function [ArrayOut] = CalcNGLDM(NGLDM,ROI)
% -------------------------------------------------------------------------
% function [ArrayOut] = CalcNGLDM(NGLDM,ROI)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates 17 NGLDM features for both 2D and 3D. 
% -------------------------------------------------------------------------
% INPUTS:
% - NGLDM: Grey level cooccurance matrix
% - ROI: the ROI from which the NGLDM is calcualted. 
% -------------------------------------------------------------------------
% OUTPUTS:
% - ArrayOut: an array of 17 GLCM fetures
% -------------------------------------------------------------------------
% AUTHOR(S):    
% - Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2015
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

nG = size(NGLDM,1);             % Number of Gray Levels
nN = size(NGLDM,2);             % Number of dependence + 1 (j = k+1)
nS = sum(NGLDM(:));             % The number of neighbourhoods.
Si = sum(NGLDM , 2);            % The number of neighbourhoods with discretised grey level i, essentially constituting a grey level histogram of the volume or slice
Sj = sum(NGLDM , 1);            % The number of neighbourhoods with dependence j, regardless of grey level.
nV = numel(find(~isnan(ROI)));  % Number of all voxels in the ROI
% FeatMatrixout = [];
% for k = 1:size(NGLDM,3)
%     tmpNGLDM = squeeze(NGLDM(:,:,k));
%     nS = sum(tmpNGLDM(:));             % Sum of all elements in NGLDM
% %     nV = numel(find(~isnan(ROI(:,:,k))));  % Number of all voxels in the ROI
%     FeatMatrixout = cat(2,FeatMatrixout , NGLDMfeatHandler(tmpNGLDM,nG,nN,nS,nV)');
% end  
% end
% 
% 
% function [ArrayOut] = NGLDMfeatHandler(NGLDM,nG,nN,nS,nV)

% Begin calling functions
ArrayOut(1)     = LowDepEmph(Sj,nN,nS);
ArrayOut(2)     = HighDepEmph(Sj,nN,nS);
ArrayOut(3)     = LowGLCountRmph(Si,nG,nS);
ArrayOut(4)     = HighGLCountEmph(Si,nG,nS);
ArrayOut(5)     = LowDepLowGLEmph(NGLDM,nG,nN,nS);
ArrayOut(6)     = LowDepHighRL(NGLDM,nG,nN,nS);
ArrayOut(7)     = HighDepLowGLEmph(NGLDM,nG,nN,nS);
ArrayOut(8)     = HighDepHighGLEmph(NGLDM,nG,nN,nS);
ArrayOut(9)     = GLnonUnif(Si,nS);
ArrayOut(10)    = GLnonUnifNormzd(Si,nS);
ArrayOut(11)    = DepCountNonUnif(Sj,nS);
ArrayOut(12)    = DepCountNonUnifNormzd(Sj,nS);
ArrayOut(13)    = DepCountPercentage(nS,nV);
ArrayOut(14)    = GLVar(NGLDM,nG,nN,nS);
ArrayOut(15)    = DepCountVar(NGLDM,nG,nN,nS);
ArrayOut(16)    = DepCountEntropy(NGLDM,nS);
ArrayOut(17)    = DepCountEnergy(NGLDM,nS);




end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature Calculation functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (LDE) Low dependence emphasis (1)
function [f_LDE] = LowDepEmph(Sj,nN,nS)
% This feature emphasises short run lengths

J = 1:nN;
tmp = Sj ./ J.^2;
f_LDE = sum(tmp)/nS;
end

%% (HDE) High dependence emphasis (2)
function [f_HDE] = HighDepEmph(Sj,nN,nS)
% This feature emphasises long run lengths

J = 1:nN;
tmp = Sj .* J.^2;
f_HDE = sum(tmp)/nS;
end

%% (LGLCE) Low grey level count emphasis (3)
function [f_LGLRE] = LowGLCountRmph(Si,nG,nS)
% This feature is a grey level analogue to f_SRE (Chu et al., 1990).
% Instead of low run lengths, low grey levels are emphasised

I = (1:nG)';
tmp = Si ./ I.^2;
f_LGLRE = sum(tmp)/nS;
end

%% (HGLCE) High grey level count emphasis (4)
function [f_HGLRE] = HighGLCountEmph(Si,nG,nS)
% The high grey level run emphasis feature is a grey level analogue to
% Frlm:lre (Chu et al., 1990). The feature emphasises high grey levels.

I = (1:nG)';
tmp = Si .* I.^2;
f_HGLRE = sum(tmp)/nS;
end

%% (LDLGLE) Low dependence low grey level emphasis (5)
function [f_SRLGLE] = LowDepLowGLEmph(NGLDM,nG,nN,nS)
% This feature emphasises runs in the upper left quadrant of the NGLDM,
% where short run lengths and low grey levels are located

[I , J] = ndgrid(1:nG , 1:nN);
tmp = NGLDM ./ (I.^2 .* J.^2);
f_SRLGLE = sum(tmp(:))/nS;
end

%% (LDHGLE) Low dependence high grey level emphasis (6)
function [f_SRHGLE] = LowDepHighRL(NGLDM,nG,nN,nS)
% This feature emphasises runs in the lower left quadrant of the NGLDM,
% where short run lengths and high grey levels are located

[I , J] = ndgrid(1:nG , 1:nN);
tmp = I.^2 .* NGLDM ./ J.^2;
f_SRHGLE = sum(tmp(:))/nS;
end

%% (HDLGLE) High dependence low grey level emphasis (7)
function [f_LRLGLE] = HighDepLowGLEmph(NGLDM,nG,nN,nS)
% This feature emphasises runs in the upper right quadrant of the NGLDM,
% where long run lengths and low grey levels are located 

[I , J] = ndgrid(1:nG , 1:nN);
tmp = J.^2 .* NGLDM ./ I.^2;
f_LRLGLE = sum(tmp(:))/nS;
end

%% (HDHGLE) High dependence high grey level emphasis (8)
function [f_LRHGLE] = HighDepHighGLEmph(NGLDM,nG,nN,nS)
% This feature emphasises runs in the lower right quadrant of the NGLDM,
% where long run lengths and high grey levels are located  

[I , J] = ndgrid(1:nG , 1:nN);
tmp = (J.^2) .* (I.^2) .* NGLDM;
f_LRHGLE = sum(tmp(:))/nS;
end

%% (GLNU) Grey level non-uniformity (9)
function [f_GLNU] = GLnonUnif(Si,nS)
% This feature assesses the distribution of runs over the grey values
% (Galloway, 1975). The feature value is low when runs are equally
% distributed along grey levels. 

f_GLNU = sum(Si.^2) / nS;
end

%% (GLNUN) Grey level non-uniformity normalized (10)
function [f_GLNUN] = GLnonUnifNormzd(Si,nS)
% This is a normalised version of the grey level non-uniformity feature.

f_GLNUN = sum(Si.^2) / nS.^2;
end

%% (RLNU) Dependence count non-uniformity (11)
function [f_RLNU] = DepCountNonUnif(Sj,nS)
% This features assesses the distribution of runs over the run lengths
% (Galloway, 1975). The feature value is low when runs are equally
% distributed along run lengths. 

f_RLNU = sum(Sj.^2) / nS;
end

%% (RLNUN) Dependence count non-uniformity normalised (12)
function [f_RLNUN] = DepCountNonUnifNormzd(Sj,nS)
% This is normalised version of the run length non-uniformity feature.

f_RLNUN = sum(Sj.^2) / nS.^2;
end

%% (DCP) Dependence count percentage (13)
function [f_RP] = DepCountPercentage(nS,nV)
% This feature assesses the fraction of the number of realised runs and the
% maximum number of potential runs (Galloway, 1975). Strongly linear or
% highly uniform ROI volumes produce a low run percentage.  

f_RP = nS / nV;
end

%% (GLV) Grey level variance (14)
function [f_GLV] = GLVar(NGLDM,nG,nN,nS)
% This feature estimates the variance in runs for the grey levels.

Pij = NGLDM / nS;
[I , ~] = ndgrid(1:nG , 1:nN);
mu = sum(sum(I .* Pij));
tmp = (I - mu).^2 .* Pij;
f_GLV = sum(tmp(:));
end

%% (DCV) Dependence count variance (15)
function [f_RLV] = DepCountVar(NGLDM,nG,nN,nS)
% This feature estimates the variance in runs for run lengths.

Pij = NGLDM / nS;
[~ , J] = ndgrid(1:nG , 1:nN);
mu = sum(sum(J .* Pij));
tmp = (J - mu).^2 .* Pij;
f_RLV = sum(tmp(:));
end

%% (DCE) Dependence count entropy (16)
function [f_RE] = DepCountEntropy(NGLDM,nS)
% It calculates the run entropy.

Pij = NGLDM / nS;
tmp = Pij .* log2(Pij + realmin);
f_RE = -sum(tmp(:));
end

%% (DCEnrg) Dependence count Energy (16)
function [f_Enrg] = DepCountEnergy(NGLDM,nS)
% It calculates the run entropy.

Pij = NGLDM / nS;
f_Enrg = sum(Pij(:).^2);
end



