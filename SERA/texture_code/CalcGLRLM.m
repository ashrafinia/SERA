function [FeatMatrixout] = CalcGLRLM(GLRLM,ROI)
% -------------------------------------------------------------------------
% function [SM_f, SS_f] = getGLCM2Dtex(ROIonly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates 16 GLRLM features for both 2D and 3D. 
% -------------------------------------------------------------------------
% INPUTS:
% - GLRLM: Grey level run length matrix
% -------------------------------------------------------------------------
% OUTPUTS:
% - FeatMatrixout: an array of 16 GLRLM fetures
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

nG = size(GLRLM,1);             % Number of Gray Levels
nR = size(GLRLM,2);             % Number of runs
% nS = sum(GLRLM(:));             % Sum of all elements in GLRLM
nV = numel(find(~isnan(ROI)));  % Number of all voxels in the ROI
FeatMatrixout = [];
for k = 1:size(GLRLM,3)
    tmpGLRLM = squeeze(GLRLM(:,:,k));
    nS = sum(tmpGLRLM(:));             % Sum of all elements in GLRLM
%     nV = numel(find(~isnan(ROI(:,:,k))));  % Number of all voxels in the ROI
    FeatMatrixout = cat(2,FeatMatrixout , GLRLMfeatHandler(tmpGLRLM,nG,nR,nS,nV)');
end  
end


function [ArrayOut] = GLRLMfeatHandler(GLRLM,nG,nR,nS,nV)

% Begin calling functions
ArrayOut(1)     = ShortRunEmph(GLRLM,nG,nR,nS);
ArrayOut(2)     = LongRunEmph(GLRLM,nG,nR,nS);
ArrayOut(3)     = LowGLRunRmph(GLRLM,nG,nR,nS);
ArrayOut(4)     = HighGLRunEmph(GLRLM,nG,nR,nS);
ArrayOut(5)     = ShortRunLowGLEmph(GLRLM,nG,nR,nS);
ArrayOut(6)     = ShortRunHighRL(GLRLM,nG,nR,nS);
ArrayOut(7)     = LongRunLowGLEmph(GLRLM,nG,nR,nS);
ArrayOut(8)     = LongRunHighGLEmph(GLRLM,nG,nR,nS);
ArrayOut(9)     = GLnonUnif(GLRLM,nS);
ArrayOut(10)    = GLnonUnifNormzd(GLRLM,nS);
ArrayOut(11)    = RunLengthNonUnif(GLRLM,nS);
ArrayOut(12)    = RunLengthNonUnifNormzd(GLRLM,nS);
ArrayOut(13)    = RunPercentage(nS,nV);
ArrayOut(14)    = GLVar(GLRLM,nG,nR,nS);
ArrayOut(15)    = RunLengthVar(GLRLM,nG,nR,nS);
ArrayOut(16)    = RunEntropy(GLRLM,nS);




end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature Calculation functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Short runs emphasis (1)
function [f_SRE] = ShortRunEmph(GLRLM,nG,nR,nS)
% This feature emphasises short run lengths

Rj = sum(GLRLM , 1);
J = 1:nR;
tmp = Rj ./ J.^2;
f_SRE = sum(tmp)/nS;
end

%% (LRE) Long runs emphasis (2)
function [f_LRE] = LongRunEmph(GLRLM,nG,nR,nS)
% This feature emphasises long run lengths

Rj = sum(GLRLM , 1);
J = 1:nR;
tmp = Rj .* J.^2;
f_LRE = sum(tmp)/nS;
end

%% (LGLRE) Low grey level run emphasis (3)
function [f_LGLRE] = LowGLRunRmph(GLRLM,nG,nR,nS)
% This feature is a grey level analogue to f_SRE (Chu et al., 1990).
% Instead of low run lengths, low grey levels are emphasised

Ri = sum(GLRLM , 2);
I = (1:nG)';
tmp = Ri ./ I.^2;
f_LGLRE = sum(tmp)/nS;
end

%% (HGLRE) High grey level run emphasis (4)
function [f_HGLRE] = HighGLRunEmph(GLRLM,nG,nR,nS)
% The high grey level run emphasis feature is a grey level analogue to
% Frlm:lre (Chu et al., 1990). The feature emphasises high grey levels.

Ri = sum(GLRLM , 2);
I = (1:nG)';
tmp = Ri .* I.^2;
f_HGLRE = sum(tmp)/nS;
end

%% (SRLGLE) Short run low grey level emphasis (5)
function [f_SRLGLE] = ShortRunLowGLEmph(GLRLM,nG,nR,nS)
% This feature emphasises runs in the upper left quadrant of the GLRLM,
% where short run lengths and low grey levels are located

[I , J] = ndgrid(1:nG , 1:nR);
tmp = GLRLM ./ (I.^2 .* J.^2);
f_SRLGLE = sum(tmp(:))/nS;
end

%% (SRHGLE) Short run high grey level emphasis (6)
function [f_SRHGLE] = ShortRunHighRL(GLRLM,nG,nR,nS)
% This feature emphasises runs in the lower left quadrant of the GLRLM,
% where short run lengths and high grey levels are located

[I , J] = ndgrid(1:nG , 1:nR);
tmp = I.^2 .* GLRLM ./ J.^2;
f_SRHGLE = sum(tmp(:))/nS;
end

%% (LRLGLE) Long run low grey level emphasis (7)
function [f_LRLGLE] = LongRunLowGLEmph(GLRLM,nG,nR,nS)
% This feature emphasises runs in the upper right quadrant of the GLRLM,
% where long run lengths and low grey levels are located 

[I , J] = ndgrid(1:nG , 1:nR);
tmp = J.^2 .* GLRLM ./ I.^2;
f_LRLGLE = sum(tmp(:))/nS;
end

%% (LRHGLE)Long run high grey level emphasis (8)
function [f_LRHGLE] = LongRunHighGLEmph(GLRLM,nG,nR,nS)
% This feature emphasises runs in the lower right quadrant of the GLRLM,
% where long run lengths and high grey levels are located  

[I , J] = ndgrid(1:nG , 1:nR);
tmp = (J.^2) .* (I.^2) .* GLRLM;
f_LRHGLE = sum(tmp(:))/nS;
end

%% (GLNU) Grey level non-uniformity (9)
function [f_GLNU] = GLnonUnif(GLRLM,nS)
% This feature assesses the distribution of runs over the grey values
% (Galloway, 1975). The feature value is low when runs are equally
% distributed along grey levels. 

Ri = sum(GLRLM , 2);
f_GLNU = sum(Ri.^2) / nS;
end

%% (GLNUN) Grey level non-uniformity normalized (10)
function [f_GLNUN] = GLnonUnifNormzd(GLRLM,nS)
% This is a normalised version of the grey level non-uniformity feature.

Ri = sum(GLRLM , 2);
f_GLNUN = sum(Ri.^2) / nS.^2;
end

%% (RLNU) Run length non-uniformity (11)
function [f_RLNU] = RunLengthNonUnif(GLRLM,nS)
% This features assesses the distribution of runs over the run lengths
% (Galloway, 1975). The feature value is low when runs are equally
% distributed along run lengths. 

Ri = sum(GLRLM , 1);
f_RLNU = sum(Ri.^2) / nS;
end

%% (RLNUN) Run length non-uniformity normalised (12)
function [f_RLNUN] = RunLengthNonUnifNormzd(GLRLM,nS)
% This is normalised version of the run length non-uniformity feature.

Ri = sum(GLRLM , 1);
f_RLNUN = sum(Ri.^2) / nS.^2;
end

%% (RP) Run percentage (13)
function [f_RP] = RunPercentage(nS,nV)
% This feature assesses the fraction of the number of realised runs and the
% maximum number of potential runs (Galloway, 1975). Strongly linear or
% highly uniform ROI volumes produce a low run percentage.  

f_RP = nS / nV;
end

%% (GLV) Grey level variance (14)
function [f_GLV] = GLVar(GLRLM,nG,nR,nS)
% This feature estimates the variance in runs for the grey levels.

Pij = GLRLM / nS;
[I , ~] = ndgrid(1:nG , 1:nR);
mu = sum(sum(I .* Pij));
tmp = (I - mu).^2 .* Pij;
f_GLV = sum(tmp(:));
end

%% (RLV) Run length variance (15)
function [f_RLV] = RunLengthVar(GLRLM,nG,nR,nS)
% This feature estimates the variance in runs for run lengths.

Pij = GLRLM / nS;
[~ , J] = ndgrid(1:nG , 1:nR);
mu = sum(sum(J .* Pij));
tmp = (J - mu).^2 .* Pij;
f_RLV = sum(tmp(:));
end

%% (RE) Run entropy (16)
function [f_RE] = RunEntropy(GLRLM,nS)
% It calculates the run entropy.

Pij = GLRLM / nS;
tmp = Pij .* log2(Pij + realmin);
f_RE = -sum(tmp(:));
end



