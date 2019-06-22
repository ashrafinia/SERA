function [FeatMatrixout] = CalcGLCM(GLCM)
% -------------------------------------------------------------------------
% function [SM_f, SS_f] = getGLCM2Dtex(ROIonly,levels)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function calculates 27 GLCM features for both 2D and 3D. 
% -------------------------------------------------------------------------
% INPUTS:
% - GLCM: Grey level cooccurance matrix
% -------------------------------------------------------------------------
% OUTPUTS:
% - FeatMatrixout: an array of 27 GLCM fetures
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

nG = size(GLCM,1);
FeatMatrixout = [];
for k = 1:size(GLCM,3)
    tmpGLCM = squeeze(GLCM(:,:,k));
    FeatMatrixout = cat(2,FeatMatrixout , GLCMfeatHandler(tmpGLCM,nG)');
end  
end


function [ArrayOut] = GLCMfeatHandler(GLCM,nG)

% Begin calling functions
ArrayOut(1)     = MaxProb(GLCM);
ArrayOut(2)     = JointAvg(GLCM,nG);
ArrayOut(3)     = JointVar(GLCM,nG,ArrayOut(2));
ArrayOut(4)     = Entropy(GLCM);
ArrayOut(5)     = DiffAvg(GLCM,nG);
ArrayOut(6)     = DiffVar(GLCM,nG,ArrayOut(5));
ArrayOut(7)     = DiffEnt(GLCM);
ArrayOut(8)     = SumAvg(GLCM,nG);
ArrayOut(9)     = SumVar(GLCM,nG,ArrayOut(8));
ArrayOut(10)    = SumEnt(GLCM);
ArrayOut(11)    = Energy(GLCM);
ArrayOut(12)    = Contrast(GLCM,nG);
ArrayOut(13)    = Dissimilarity(GLCM,nG);
ArrayOut(14)    = InvDiff(GLCM,nG);
ArrayOut(15)    = InvDiffNorm(GLCM,nG);
ArrayOut(16)    = InvDiffMom(GLCM,nG);
ArrayOut(17)    = InvDiffMomNorm(GLCM,nG);
ArrayOut(18)    = InvVar(GLCM,nG);
ArrayOut(19)    = Correlation(GLCM,nG);
ArrayOut(20)    = AutoCorr(GLCM,nG);
ArrayOut(21)    = ClusterTend(GLCM,nG);
ArrayOut(22)    = ClusterShade(GLCM,nG);
ArrayOut(23)    = ClusterProm(GLCM,nG);
ArrayOut(24)    = InfoCorr1(GLCM,ArrayOut(4),nG);
ArrayOut(25)    = InfoCorr2(GLCM,ArrayOut(4),nG);
% ArrayOut(26)    = Agreement(GLCM);
% ArrayOut(27)    = SumAvg2(GLCM,nG);




end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature Calculation functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Max Prob(1)
function [f_joint_max] = MaxProb(GLCM)
% The joint average is the grey level weighted sum of joint probabilities

f_joint_max = max(GLCM(:));
end

%% Joint Average (2)
function [f_joint_avg] = JointAvg(GLCM,nG)
% The joint average is the grey level weighted sum of joint probabilities

tmp = GLCM .* repmat((1:nG)' , 1, nG);
f_joint_avg = sum(tmp(:));
end

%% Joint variance (3)
function [f_joint_var] = JointVar(GLCM,nG , mu)
% The joint variance, which is also called sum of squares

tmp = GLCM .* repmat(((1:nG)' - mu).^2 , 1, nG);
f_joint_var = sum(tmp(:));
end

%% Entropy (4)
function [f_entropy] = Entropy(GLCM)
% The joint Entropy

tmp = GLCM .* log2(GLCM+realmin);
f_entropy = -sum(tmp(:));
end

%% Diference average (5)
function [f_diffavg] = DiffAvg(GLCM,nG)
% The average for the diagonal probabilities is defned as:

pDiag = DiagProb(GLCM);
tmp = (0:(nG-1)) .* pDiag;
f_diffavg = sum(tmp(:));
end

%% Diference variance (6)
function [f_diffvar] = DiffVar(GLCM,nG,muDiffAvg)
% The variance for the diagonal probabilities

pDiag = DiagProb(GLCM);
tmp = ((0:(nG-1))-muDiffAvg).^2 .* pDiag;
f_diffvar = sum(tmp(:));
end

%% Difference Entropy (7)
function [f_diffent] = DiffEnt(GLCM)
% The entropy for the diagonal probabilities

pDiag = DiagProb(GLCM);
tmp = pDiag .* log2(pDiag+realmin);
f_diffent = -sum(tmp(:));
end

%% Sum average (8)
function [f_sumavg] = SumAvg(GLCM,nG)
% The average for the diagonal probabilities is defned as:

pCross = CrossProb(GLCM);
tmp = (2:(2*nG)) .* pCross;
f_sumavg = sum(tmp(:));
end

%% Sum variance (9)
function [f_sumvar] = SumVar(GLCM,nG,muSumAvg)
% The variance for the diagonal probabilities

pCross = CrossProb(GLCM);
tmp = ((2:(2*nG))-muSumAvg).^2 .* pCross;
f_sumvar = sum(tmp(:));
end

%% Sum Entropy (10)
function [f_sument] = SumEnt(GLCM)
% The entropy for the diagonal probabilities

pCross = CrossProb(GLCM);
tmp = pCross .* log2(pCross+realmin);
f_sument = -sum(tmp(:));
end

%% Energy (11)
function [f_nrg] = Energy(GLCM)
% The angular second moment which represent the energy

f_nrg = sum(GLCM(:).^2);
end

%% Contrast (12)
function [f_cont] = Contrast(GLCM,nG)
% Contrast assesses grey level variations (Haralick et al., 1973). Hence 
% elements that represent large grey level diferences receive greater weight

[I , J ] = meshgrid(1:nG , 1:nG);
T = I-J;
tmp = (T .^ 2) .* GLCM;
f_cont = sum(tmp(:));
end

%% Dissimilarity (13)
function [f_diss] = Dissimilarity(GLCM,nG)
% The entropy for the diagonal probabilities

[I , J ] = meshgrid(1:nG , 1:nG);
T = I-J;
tmp = abs(T) .* GLCM;
f_diss = sum(tmp(:));
end

%% Homogeneity (Inverse Difference) (14)
function [f_invdiff] = InvDiff(GLCM,nG)
% Inverse diference is a measure of homogeneity. Grey level co-occurrences
% with a large diference in levels are weighed less, thus lowering the
% total feature score. The feature score is maximal if all grey levels are
% the same. 

[I , J ] = meshgrid(1:nG , 1:nG);
T = I-J;
tmp = GLCM ./ (1+abs(T));
f_invdiff = sum(tmp(:));
end

%% Normalized Homogeneity (Inverse diference normalised) (15)
function [f_invdiffN] = InvDiffNorm(GLCM,nG)
% Clausi (Clausi, 2002) suggests normalising inverse diference to improve
% classifcation ability of this feature.

[I , J ] = meshgrid(1:nG , 1:nG);
T = I-J;
tmp = GLCM ./ (1+abs(T)/nG);
f_invdiffN = sum(tmp(:));
end

%% (Homogeneity 2) Inverse diference moment (16)
function [f_invdiffM] = InvDiffMom(GLCM,nG)
% Inverse diference moment (Haralick et al., 1973) is similar in concept to
% the inverse diference feature, but with lower weights for elements that
% are further from the diagonal.  

[I , J ] = meshgrid(1:nG , 1:nG);
T = I-J;
tmp = GLCM ./ (1+T.^2);
f_invdiffM = sum(tmp(:));
end

%% Inverse diference moment normalised (17)
function [f_invdiffMN] = InvDiffMomNorm(GLCM,nG)
% Clausi (Clausi, 2002) suggests normalising inverse diference moment to
% improve classifcation performance of this feature.

[I , J ] = meshgrid(1:nG , 1:nG);
T = I-J;
tmp = GLCM ./ (1+T.^2/nG^2);
f_invdiffMN = sum(tmp(:));
end

%% Inverse variance (18)
function [f_invVar] = InvVar(GLCM,nG)
% The inverse variance feature is defned as:

[I , J ] = meshgrid(1:nG , 1:nG);
T = I-J;
tmp = GLCM ./ T.^2;
tmp2 = triu(tmp,1);
f_invVar = 2* sum(tmp2(:));
end

%% Correlation (19)
function [f_corr] = Correlation(GLCM,nG)
% Calculates the correlation

Pi = sum(GLCM,2); 
Pj = sum(GLCM,1);
Ui = sum((1:nG)' .* Pi);
Uj = sum((1:nG) .* Pj);
Si = sqrt(sum(((1:nG)'-Ui).^2 .* Pi));

[I , J ] = meshgrid(1:nG , 1:nG);
tmp = (I - Ui) .* (J - Uj) .* GLCM;
f_corr = sum(tmp(:)) / Si.^2;
f_corr(isnan(f_corr)) = 0; 
end

%% AutoCorrelation (20)
function [f_autocorr] = AutoCorr(GLCM,nG)
% Calculate autocorrelation

[I , J ] = meshgrid(1:nG , 1:nG);
tmp = I .* J .* GLCM;
f_autocorr = sum(tmp(:));
end

%% Cluster tendency (21)
function [f_clstnd] = ClusterTend(GLCM,nG)
% Cluster tendency is defned as:

Pi = sum(GLCM,2); 
Ui = sum((1:nG)' .* Pi);
[I , J ] = meshgrid(1:nG , 1:nG);
tmp = (I + J - 2*Ui).^2 .* GLCM;
f_clstnd = sum(tmp(:));
end

%% Cluster shade (22)
function [f_clsshd] = ClusterShade(GLCM,nG)
% Cluster tendency is defned as:

Pi = sum(GLCM,2); 
Ui = sum((1:nG)' .* Pi);
[I , J ] = meshgrid(1:nG , 1:nG);
tmp = (I + J - 2*Ui).^3 .* GLCM;
f_clsshd = sum(tmp(:));
end

%% Cluster prominence (23)
function [f_clsprm] = ClusterProm(GLCM,nG)
% Cluster prominence is defned as:

Pi = sum(GLCM,2); 
Ui = sum((1:nG)' .* Pi);
[I , J ] = meshgrid(1:nG , 1:nG);
tmp = (I + J - 2*Ui).^4 .* GLCM;
f_clsprm = sum(tmp(:));
end

%% First measure of information correlation (24)
function [f_ic1] = InfoCorr1(GLCM,Ent,nG)
% Information theoretic correlation is estimated using two diferent measure

Pi = sum(GLCM,2); 
Pj = sum(GLCM,1);
tmp = Pi .* log2(Pi+realmin);
HX = -sum(tmp(:));
tmp = Pj .* log2(Pj+realmin);
HY = -sum(tmp(:));
tmp = GLCM .* log2(repmat(Pi,1,nG).*repmat(Pj,nG,1)+realmin);
HXY1 = -sum(tmp(:));
f_ic1 = (Ent - HXY1) / max(HX,HY);
f_ic1(isnan(f_ic1)) = 0;
end

%% First measure of information correlation (25)
function [f_ic2] = InfoCorr2(GLCM,Ent,nG)
% Information theoretic correlation is estimated using two diferent measure

Pi = sum(GLCM,2); 
Pj = sum(GLCM,1);
tmp = (repmat(Pi,1,nG).*repmat(Pj,nG,1)) .* log2(repmat(Pi,1,nG).*repmat(Pj,nG,1)+realmin);
HXY2 = -sum(tmp(:));
if Ent > HXY2
    f_ic2 = 0;
else
    f_ic2 = sqrt(1-exp(-2*(HXY2-Ent)));
end
end

%% Agreement (26)
function [f_agr] = Agreement(GLCM)
% Calculates Agreement

Po = sum(diag(GLCM));
Pe = sum(diag(GLCM*GLCM));
f_agr = (Po-Pe) / (1-Pe);
f_agr(isnan(f_agr))= 0; 
end


%% Sum average 2 (27)
function [f_sumavg2] = SumAvg2(GLCM,nG)
% The average for the diagonal probabilities is defned as:

[I , J ] = meshgrid(1:nG , 1:nG);
tmp = I .* GLCM + J .* GLCM;
f_sumavg2 = sum(tmp(:) / nG^2);
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Extra calculation functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Diagonal probabilities
function [pDiag] = DiagProb(GLCM)
% This function calculates the diagonal probabilities. 

pDiag = zeros(1,size(GLCM,1));
for k=1:size(GLCM,1)
    if k==1
        pDiag(k) = sum(diag(GLCM,k-1));
    else
        pDiag(k) = sum(diag(GLCM,k-1)) + sum(diag(GLCM,-(k-1)));
    end
end
end

%% Cross-diagonal probability
function [pCross] = CrossProb(GLCM)
% This function calculates the Cross-diagonal probability. 

nG = size(GLCM,1);
pCross = zeros(1,nG*2-1);
GF = fliplr(GLCM);
for k=1:nG*2-1
    pCross(k) = sum(diag(GF,nG - k));
end
end


%%

% StructOut.f_joint_max	=	ArrayOut(1)	;
% StructOut.f_joint_avg	=	ArrayOut(2)	;
% StructOut.f_joint_var	=	ArrayOut(3)	;
% StructOut.f_entropy     =	ArrayOut(4)	;
% StructOut.f_diffavg     =	ArrayOut(5)	;
% StructOut.f_diffvar     =	ArrayOut(6)	;
% StructOut.f_diffent     =	ArrayOut(7)	;
% StructOut.f_sumavg      =	ArrayOut(8)	;
% StructOut.f_sumvar      =	ArrayOut(9)	;
% StructOut.f_sument      =	ArrayOut(10);
% StructOut.f_nrg         =	ArrayOut(11);
% StructOut.f_cont        =	ArrayOut(12);
% StructOut.f_diss        =	ArrayOut(13);
% StructOut.f_invdiff     =	ArrayOut(14);
% StructOut.f_invdiffN	=	ArrayOut(15);
% StructOut.f_invdiffM	=	ArrayOut(16);
% StructOut.f_invdiffMN	=	ArrayOut(17);
% StructOut.f_invVar      =	ArrayOut(18);
% StructOut.f_corr        =	ArrayOut(19);
% StructOut.f_autocorr	=	ArrayOut(20);
% StructOut.f_clstnd      =	ArrayOut(21);
% StructOut.f_clsshd      =	ArrayOut(22);
% StructOut.f_clsprm      =	ArrayOut(23);
% StructOut.f_ic1         =	ArrayOut(24);
% StructOut.f_ic2         =	ArrayOut(25);
% StructOut.f_agr         =	ArrayOut(26);
% StructOut.f_sumavg2     =	ArrayOut(27);

% OUT.f_joint_max   = MaxProb(GLCM);
% OUT.f_joint_avg   = JointAvg(GLCM,nG);
% OUT.f_joint_var   = JointVar(GLCM,nG , OUT.f_joint_avg);
% OUT.f_entropy     = Entropy(GLCM);
% OUT.f_diffavg     = DiffAvg(GLCM,nG);
% OUT.f_diffvar     = DiffVar(GLCM,nG,OUT.f_joint_avg);
% OUT.f_diffent     = DiffEnt(GLCM);
% OUT.f_sumavg      = SumAvg(GLCM,nG);
% OUT.f_sumvar      = SumVar(GLCM,nG,OUT.f_joint_avg);
% OUT.f_sument      = SumEnt(GLCM);
% OUT.f_nrg         = Energy(GLCM);
% OUT.f_cont        = Contrast(GLCM,nG);
% OUT.f_diss        = Dissimilarity(GLCM,nG);
% OUT.f_invdiff     = InvDiff(GLCM,nG);
% OUT.f_invdiffN    = InvDiffNorm(GLCM,nG);
% OUT.f_invdiffM    = InvDiffMom(GLCM,nG);
% OUT.f_invdiffMN   = InvDiffMomNorm(GLCM,nG);
% OUT.f_invVar      = InvVar(GLCM,nG);
% OUT.f_corr        = Correlation(GLCM,nG);
% OUT.f_autocorr    = AutoCorr(GLCM,nG);
% OUT.f_clstnd      = ClusterTend(GLCM,nG);
% OUT.f_clsshd      = ClusterShade(GLCM,nG);
% OUT.f_clsprm      = ClusterProm(GLCM,nG);
% OUT.f_ic1         = InfoCorr1(GLCM,OUT.f_entropy,nG);
% OUT.f_ic2         = InfoCorr2(GLCM,OUT.f_entropy,nG);
% OUT.f_agr         = Agreement(GLCM);
% OUT.f_sumavg2     = SumAvg2(GLCM,nG);
