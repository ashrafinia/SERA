function [textures] = getNGTDMtextures(NGTDM,countValid,Aarray)
% -------------------------------------------------------------------------
% function [textures] = getNGTDMtextures(NGTDM,countValid)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes texture features from an input Neighborhood
% Gray-Tone Difference Matrix (NGTDM).
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Amadasun, M., & King, R. (1989). Textural Features Corresponding to 
%     Textural Properties. IEEE Transactions on Systems Man and Cybernetics,
%     19(5), 1264–1274.
% -------------------------------------------------------------------------
% INPUTS:
% - NGTDM: Neighborhood Gray-Tone Difference Matrix.
% - countValid: Number of valid voxels used in the NGTDM computation. 
%               Required for the computation of texture features in 
%               'getNGTDMtextures.m'
%
% ** 'NGTDM' and 'countValid' should be outputs from 'getNGTDM.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different NGTDM texture
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


% PRELIMINARY
nV = sum(countValid);
Pi = countValid./nV; % Now representing the probability of gray-level occurences
nG = length(NGTDM);
nGp = sum(Pi~=0);
pValid = find(Pi>0);
nValid = length(pValid);
% Si = accumarray(Aarray(:,1),abs(Aarray(:,1) - Aarray(:,2)));


% COMPUTATION OF TEXTURES
% 1. Coarseness, Ref.[1]
textures.Coarseness = (((Pi')*NGTDM) + eps)^(-1);
textures.Coarseness = min(textures.Coarseness , 10^6);

% 2. Contrast, Ref.[1]
val = 0;
for i = 1:nG
    for j = 1:nG
        val = val + Pi(i)*Pi(j)*(i-j)^2;
    end
end
textures.Contrast = val*sum(NGTDM)/(nGp*(nGp-1)*nV);

% 3. Busyness, Ref.[1]
denom = 0;
for i = 1:nValid
    for j = 1:nValid
        denom = denom + abs(pValid(i)*Pi(pValid(i))-pValid(j)*Pi(pValid(j)));
    end
end
textures.Busyness = ((Pi')*NGTDM)/denom;
if isinf(textures.Busyness), textures.Busyness=0; end

% 4. Complexity, Ref.[1]
val = 0;
for i = 1:nValid
    for j = 1:nValid
        val = val + (abs(pValid(i)-pValid(j))/(nV*(Pi(pValid(i)) + Pi(pValid(j)))))*(Pi(pValid(i))*NGTDM(pValid(i)) + Pi(pValid(j))*NGTDM(pValid(j)));
    end
end
textures.Complexity = val;

% 5. Strength, Ref.[1]
val = 0;
for i = 1:nValid
    for j = 1:nValid
        val = val + (Pi(pValid(i))+Pi(pValid(j)))*(pValid(i)-pValid(j))^2;
    end
end
textures.Strength = val/(eps+sum(NGTDM));

end