function [textures] = getGLCMtextures(GLCM)
% -------------------------------------------------------------------------
% function [textures] = getGLCMtextures(GLCM))
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes texture features from an input Gray-Level 
% Co-occurence Matrix (GLCM).
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Haralick, R. M., Shanmugam, K., & Dinstein, I. (1973). Textural 
%     features for image classification. IEEE Transactions on Systems, Man 
%     and Cybernetics, smc 3(6), 610–621.
% [2] Assefa, D., Keller, H., Ménard, C., Laperriere, N., Ferrari, R. J., & 
%     Yeung, I. (2010). Robust texture features for response monitoring of 
%     glioblastoma multiforme on T1 -weighted and T2 -FLAIR MR images: A 
%     preliminary investigation in terms of identification and segmentation. 
%     Medical Physics, 37(4), 1722–1736.
% [3] Thibault, G. (2009). Indices de formes et de textures: de la 2D vers 
%     la 3D. Application au classement de noyaux de cellules. PhD Thesis, 
%     Université AIX-Marseille: p.172.
% -------------------------------------------------------------------------
% INPUTS:
% - GLCM: Gray-Level Co-occurence Matrix.
%
% ** 'GLCM' should be the output from 'getGLCM.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different GLCM texture
%             features as defined below.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres, <mart.vallieres@gmail.com>
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
textures = struct;
matrixtemp = GLCM;
GLCM = GLCM/(sum(GLCM(:))); % Normalization of GLCM
nL = max(size(GLCM));
indVect = 1:nL;
[colGrid,rowGrid] = meshgrid(indVect,indVect);


% COMPUTATION OF TEXTURE FEATURES
% 1. Energy, Ref.[1]
textures.Energy = sum(sum(GLCM.^2));

% 2. Contrast, Ref.[1]
contrast = 0.0;
for n = 0:nL-1%sum start from 1 based eq contrast
   temp = 0;
   for i = 1:nL
      for j = 1:nL
         if (abs(i-j) == n)
            temp = temp+GLCM(i,j);
         end
      end
   end
   contrast = contrast + n^2*temp;
end
textures.Contrast = contrast;

% 3. Entropy, Ref.[1]
textures.Entropy = -sum(sum(GLCM.*log2(GLCM + realmin)));

% 4. Homogeneity, adapted from Ref.[1]
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + GLCM(i,j)/(1+abs(i-j));
   end
end
textures.Homogeneity = temp;

temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + GLCM(i,j)/(1+abs(i-j).^2);
   end
end
textures.Homogeneity2 = temp;

% 5. Correlation, adapted from Ref. [1] (this definition from MATLAB is preferred from the original one in [1])
textures.Correlation = graycoprops(round(matrixtemp),'Correlation');
textures.Correlation = struct2cell(textures.Correlation);
textures.Correlation = textures.Correlation{1};

%new correlation
px = sum(GLCM,2);% one row add as a number
py = sum(GLCM);
ux = mean(px);
uy = mean(py);
sigx = std(px);
sigy = std(py);
temp =0;
for i = 1:nL
    for j = 1:nL
        temp = temp + i*j*GLCM(i,j);        
    end
end
textures.Correlation2 = (temp - ux*uy)/(sigx*sigy);


% new1  autocorrelation
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + i*j*GLCM(i,j);
   end
end
textures.Autocorrelation = temp;

%new2 Cluster prominence
ux = mean(GLCM,2);
uy = mean(GLCM);
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + (i+j-ux(i)-uy(j)).^4.*GLCM(i,j);
   end
end
textures.Clusterprominence = temp;

%new3 CLuster shade
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + (i+j-ux(i)-uy(j)).^3.*GLCM(i,j);
   end
end
textures.Clustershade = temp;

%new4 Cluster tendency
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + (i+j-ux(i)-uy(j)).^2.*GLCM(i,j);
   end
end
textures.Clustertendency = temp;

%new5 Difference entropy
temp = zeros(1,nL);
for k = 0:nL-1
  for i = 1:nL
     for j = 1:nL
       if(i-j==k)
        temp(k+1) = temp(k+1) + GLCM(i,j);%noly i>=j those elements were sumed
        else
        end
     end
  end
end
px_y = temp;
textures.Differenceentropy = -px_y*(log2(px_y+realmin))';

% sum entropy
temp = zeros(1,2*nL-2+1);
for k = 2:2*nL
  for i = 1:nL
     for j = 1:nL
       if(i+j==k)
        temp(k-1) = temp(k-1) + GLCM(i,j);%noly i>=j those elements were sumed
        else
        end
     end
  end
end
pxy = temp;
textures.Sumentropy = -pxy*(log2(pxy+realmin))';

%new sum average
i = 2:2*nL;
textures.SumAverage2 = i*pxy';

%new variance
u = mean(GLCM(:));
temp = 0;
for i = 1:nL;
    for j = 1:nL;
        temp = temp + (i-u)^2*GLCM(i,j);
    end
end
textures.Variance2 = temp;

%sum variance
SE = textures.Sumentropy;
temp = 0;
for i = 1:length(pxy)
    temp = temp + (i+1-SE).^2*pxy(i);
end
textures.Sumvariance = temp;

%new6 informational measure of correlation1(ICM1and ICM2)
H = textures.Entropy;
px = sum(GLCM,2);%sum of the row
py = sum(GLCM);%sum of the clumn
HXY1 = 0;
HXY2 = 0;
for i = 1:nL
   for j = 1:nL
      HXY1 = HXY1-GLCM(i,j)*log2(px(i)*py(j)+realmin);
      HXY2 = HXY2-px(i)*py(j)*log2(px(i)*py(j)+realmin);
   end
end
HX = -px'*(log2(px+realmin));%entropies of px and py
HY = -py*(log2(py+realmin))';
textures.ICM1 = (H-HXY1)./(max(HX,HY)+realmin);
textures.ICM2 = (1-exp(-2*(HXY2-H))).^0.5;

% %maximal correlation coefficient
% Q = zeros(nL,nL);
% for i = 1:nL;
%     for j = 1:nL;
%         for k = 1:nL
%             Q(i,j) = Q(i,j) + GLCM(i,k)*GLCM(j,k)/(px(i)*py(k));
%         end 
%     end 
% end 
% Q = Q((~isnan(Q))&(~isinf(Q)));
% Q = reshape(Q,nL,nL);
% [V,D] = eig(Q);
% D = diag(D);
% [D,I] = sort(D,'descend');
% textures.Maxcorrcoeff = D(2).^0.5;


% 6. Variance, Ref.[2]; and 7. SumAverage, Ref.[2]. (adapted from Variance and SumAverage metrics defined by Haralick in Ref. [1])
% However, in order to compare GLCMs of different sizes, the metrics
% are divided by the total number of elements in the GLCM (nL*nL). Also,
% there is probably an error in Assefa's paper [2]: in the variance equation,
% 'u' should appropriately be replaced by 'ux' and 'uy' as calculated in A1
% and A2 of the same paper (read 'ui' and 'uj' in our code).
ui = indVect*sum(GLCM,2);
uj = indVect*sum(GLCM)';
tempS = rowGrid.*GLCM + colGrid.*GLCM;
tempV = (rowGrid-ui).^2.*GLCM + (colGrid-uj).^2.*GLCM;
textures.SumAverage = 0.5*sum(tempS(:))/(nL^2);
textures.Variance = 0.5*sum(tempV(:))/(nL^2);

% 8. Dissimilarity, Ref.[3] 
diffMat = abs(rowGrid-colGrid);
temp = diffMat.*GLCM;
textures.Dissimilarity = sum(temp(:));

% 9. max possibility
textures.MaxPossibility = max(max(GLCM));

% 10. inverse difference moment
tempIDM = 0;
tempIDMN = 0;
tempIDN = 0;
for i = 1:nL
   for j = 1:nL
       tempIDMN = tempIDMN + GLCM(i,j)/(1+(i-j)^2/(nL^2));
       tempIDN = tempIDN + GLCM(i,j)/(1+(i-j)/nL);
       if i~=j
          tempIDM = tempIDM + GLCM(i,j)/((i-j)^2);
       else
       end
   end
end
textures.InvDifMoment = tempIDM;
textures.IDMN = tempIDMN;
textures.IDN = tempIDN;



% agreement
pe = 0;
p0 = trace(GLCM);
k = 0;
for i = 1:nL;
    pe = pe + GLCM(i,:)*GLCM(:,i);
end
k = (p0 - pe)/(1-pe+realmin);
textures.Agreement = k;


end








