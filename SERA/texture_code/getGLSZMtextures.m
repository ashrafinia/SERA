function [textures] = getGLSZMtextures(GLSZM)
% -------------------------------------------------------------------------
% function [textures] = getGLSZMtextures(GLSZM)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes texture features from an input Gray-Level Size
% Zone Matrix (GLSZM).
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Galloway, M. M. (1975). Texture analysis using gray level run lengths. 
%     Computer Graphics and Image Processing, 4(2), 172–179.
% [2] Chu, A., Sehgal, C. M., & Greenleaf, J. F. (1990). Use of gray value 
%     distribution of run lengths for texture analysis. Pattern Recognition
%     Letters, 11(6), 415-419.
% [3] Dasarathy, B. V., & Holder, E. B. (1991). Image characterizations 
%     based on joint gray level-run length distributions. Pattern 
%     Recognition Letters, 12(8), 497-502.
% [4] Thibault, G., Fertil, B., Navarro, C., Pereira, S., Cau, P., Levy, 
%     N., Mari, J.-L. (2009). Texture Indexes and Gray Level Size Zone 
%     Matrix. Application to Cell Nuclei Classification. In Pattern 
%     Recognition and Information Processing (PRIP) (pp. 140–145).
% -------------------------------------------------------------------------
% INPUTS:
% - GLSZM: Gray-Level Size Zone Matrix.
%
% ** 'GLSZM' should be the output from 'getGLSZM.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different GLSZM texture
%             features as defined below.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres
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


% USEFUL MATRICES, VECTORS AND QUANTITIES
sz = size(GLSZM); % Size of GLSZM
nRuns = sum(GLSZM(:));
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLSZM
pg = sum(GLSZM,2)'; % Gray-Level Run-Number Vector
pr = sum(GLSZM,1); % Run-Length Run-Number Vector


% COMPUTATION OF TEXTURE FEATURES
% 1. Small Zone Emphasis (SZE), Ref.[1,4]
textures.SZE = (pr*(cVect.^(-2))')/nRuns;

% 2. Large Zone Emphasis (LZE), Ref.[1,4]
textures.LZE = (pr*(cVect.^2)')/nRuns;


% 3. Low Gray-Level Zone Emphasis (LGZE), Ref.[2,4]
textures.LGZE = (pg*(rVect.^(-2))')/nRuns;

% 4. High Gray-Level Zone Emphasis (HGZE), Ref.[2,4]
textures.HGZE = (pg*(rVect.^2)')/nRuns;

% 5. Small Zone Low Gray-Level Emphasis (SZLGE), Ref.[3,4]
textures.SZLGE = sum(sum(GLSZM.*(rMat.^(-2)).*(cMat.^(-2))))/nRuns;

% 6. Small Zone High Gray-Level Emphasis (SZHGE), Ref.[3,4]
textures.SZHGE = sum(sum(GLSZM.*(rMat.^2).*(cMat.^(-2))))/nRuns;

% 7. Large Zone Low Gray-Level Emphasis (LZLGE), Ref.[3,4]
textures.LZLGE = sum(sum(GLSZM.*(rMat.^(-2)).*(cMat.^2)))/nRuns;

% 8. Large Zone High Gray-Level Emphasis (LZHGE), Ref.[3,4]
textures.LZHGE = sum(sum(GLSZM.*(rMat.^2).*(cMat.^2)))/nRuns;


% 9. Gray-Level Nonuniformity (GLN), adapted from Ref.[1,4]
textures.GLN = sum(pg.^2)/nRuns;

% 10. Gray-Level Nonuniformity Normalized (GLN), adapted from Ref.[1,4]
textures.GLNN = sum(pg.^2)/nRuns.^2;

% 11. Zone-Size Nonuniformity (ZSN), adapted from Ref.[1,4]
textures.ZSN = sum(pr.^2)/nRuns;

% 12. Zone-Size Nonuniformity Normalized (ZSNN), adapted from Ref.[1,4]
textures.ZSNN = sum(pr.^2)/nRuns.^2;


% 13. Zone Percentage (ZP), adapted from Ref.[1,4]
textures.ZP = nRuns/(pr*cVect');

% 14. Gray-Level Variance (GLV), adapted from Ref.[4]
mu = sum(rVect *GLSZM)/nRuns;
textures.GLV = sum((rVect-mu).^2 *GLSZM)/nRuns;

% 15. Zone-Size Variance (ZSV), adapted from Ref.[4]
mu = sum(GLSZM*cVect')/nRuns;
textures.ZSV = sum(GLSZM * (cVect'-mu).^2 )/nRuns;

% 16. Zone-size Entropy 
textures.Entropy = -sum((GLSZM(:)/nRuns) .* log2(GLSZM(:)/nRuns+realmin));

end