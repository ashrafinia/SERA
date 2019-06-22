function [out] = GLCMFeatures(glcm)
% 
% Features computed 
% Autocorrelation: [2]   
% Cluster Prominence: [2]                   
% Cluster Shade: [2] 
% Contrast: [1]                                         
% Correlation: [1]                        
% Difference entropy [1] 
% Difference variance [1]                   
% Dissimilarity: [2]                        
% Energy: [1]                    
% Entropy: [2]       
% Homogeneity: (Inverse Difference Moment) [2,1] 
% Information measure of correlation1 [1]   
% Informaiton measure of correlation2 [1]  
% Inverse difference (Homogeneity in matlab): [3]                              
% Maximum probability: [2]                    
% Sum average [1]   
% Sum entropy [1]  
% Sum of sqaures: Variance [1]    
% Sum variance [1]   
%
% References:
% 1. R. M. Haralick, K. Shanmugam, and I. Dinstein, Textural Features of
% Image Classification, IEEE Transactions on Systems, Man and Cybernetics,
% vol. SMC-3, no. 6, Nov. 1973
% 2. L. Soh and C. Tsatsoulis, Texture Analysis of SAR Sea Ice Imagery
% Using Gray Level Co-Occurrence Matrices, IEEE Transactions on Geoscience
% and Remote Sensing, vol. 37, no. 2, March 1999.
% 3. D A. Clausi, An analysis of co-occurrence texture statistics as a
% function of grey level quantization, Can. J. Remote Sensing, vol. 28, no.
% 1, pp. 45-62, 2002
%
%
% Started from Avinash Uppupuri's code on Matlab file exchange. It has then
% been vectorized. Three features were not implemented correctly in that
% code, it has since then been changed. The features are: 
%   * Sum of squares: variance
%   * Difference variance
%   * Sum Variance


if ((nargin > 1) || (nargin == 0))
    error('Too many or too few input arguments')
else
    if ((size(glcm,1) <= 1) || (size(glcm,2) <= 1))
        error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcm,1) ~= size(glcm,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end  
end

% Normalize the GLCMs
glcm = bsxfun(@rdivide,glcm,sum(sum(glcm)));

% Get size of GLCM
nGrayLevels = size(glcm,1);
nglcm = size(glcm,3);

% checked 
out.autoCorrelation                     = zeros(1,nglcm); % Autocorrelation: [2] 
out.clusterProminence                   = zeros(1,nglcm); % Cluster Prominence: [2]
out.clusterShade                        = zeros(1,nglcm); % Cluster Shade: [2]
out.contrast                            = zeros(1,nglcm); % Contrast: matlab/[1,2]
out.correlation                         = zeros(1,nglcm); % Correlation: [1,2]
out.differenceEntropy                   = zeros(1,nglcm); % Difference entropy [1]
out.differenceVariance                  = zeros(1,nglcm); % Difference variance [1]
out.dissimilarity                       = zeros(1,nglcm); % Dissimilarity: [2]
out.energy                              = zeros(1,nglcm); % Energy: matlab / [1,2]
out.entropy                             = zeros(1,nglcm); % Entropy: [2]
out.homogeneity                         = zeros(1,nglcm); % Homogeneity: [2] (inverse difference moment)
out.informationMeasureOfCorrelation1    = zeros(1,nglcm); % Information measure of correlation1 [1]
out.informationMeasureOfCorrelation2    = zeros(1,nglcm); % Informaiton measure of correlation2 [1]
out.inverseDifference                   = zeros(1,nglcm); % Homogeneity in matlab
% out.inverseDifferenceMomentNormalized   = zeros(1,nglcm); % Normalized Homogeneity
% out.inverseDifferenceNormalized         = zeros(1,nglcm); % Normalized inverse difference
out.maximumProbability                  = zeros(1,nglcm); % Maximum probability: [2]
out.sumAverage                          = zeros(1,nglcm); % Sum average [1]    
out.sumEntropy                          = zeros(1,nglcm); % Sum entropy [1]
out.sumOfSquaresVariance                = zeros(1,nglcm); % Sum of sqaures: Variance [1]
out.sumVariance                         = zeros(1,nglcm); % Sum variance [1]

glcmMean = zeros(nglcm,1);
uX = zeros(nglcm,1);
uY = zeros(nglcm,1);
sX = zeros(nglcm,1);
sY = zeros(nglcm,1);

% pX pY pXplusY pXminusY
pX = zeros(nGrayLevels,nglcm); % Ng x #glcms[1]  
pY = zeros(nGrayLevels,nglcm); % Ng x #glcms[1]
pXplusY = zeros((nGrayLevels*2 - 1),nglcm); %[1]
pXminusY = zeros((nGrayLevels),nglcm); %[1]
% HXY1 HXY2 HX HY
HXY1 = zeros(nglcm,1);
HX   = zeros(nglcm,1);
HY   = zeros(nglcm,1);
HXY2 = zeros(nglcm,1);

% Create indices for vectorising code:
sub   = 1:nGrayLevels*nGrayLevels;
[I,J] = ind2sub([nGrayLevels,nGrayLevels],sub);

% Loop over all GLCMs
for k = 1:nglcm 
    currentGLCM = glcm(:,:,k);
    glcmMean(k) = mean2(currentGLCM);
    
    % For symmetric GLCMs, uX = uY
    uX(k)   = sum(I.*currentGLCM(sub));
    uY(k)   = sum(J.*currentGLCM(sub));
    sX(k)   = sum((I-uX(k)).^2.*currentGLCM(sub));
    sY(k)   = sum((J-uY(k)).^2.*currentGLCM(sub));

    out.contrast(k)             = sum(abs(I-J).^2.*currentGLCM(sub)); %OK
    out.dissimilarity(k)        = sum(abs(I - J).*currentGLCM(sub)); %OK
    out.energy(k)               = sum(currentGLCM(sub).^2); % OK
    out.entropy(k)              = -nansum(currentGLCM(sub).*log(currentGLCM(sub))); %OK
    out.inverseDifference(k)    = sum(currentGLCM(sub)./( 1 + abs(I-J) )); %OK
    out.homogeneity(k)          = sum(currentGLCM(sub)./( 1 + (I - J).^2)); %OK
    
%     out.inverseDifferenceNormalized(k)      = sum(currentGLCM(sub)./( 1 + abs(I-J)/nGrayLevels )); %OK
%     out.inverseDifferenceMomentNormalized(k)= sum(currentGLCM(sub)./( 1 + ((I - J)/nGrayLevels).^2)); %OK

    out.sumOfSquaresVariance(k) = sum(currentGLCM(sub).*((I - uX(k)).^2)); %<----- N.B! Wrong implementation previously!!
    out.maximumProbability(k)   = max(currentGLCM(:));
    
    pX(:,k) = sum(currentGLCM,2); %OK
    pY(:,k) = sum(currentGLCM,1)'; %OK
    
    tmp1 = [(I+J)' currentGLCM(sub)'];
    tmp2 = [abs((I-J))' currentGLCM(sub)'];
    idx1 = 2:2*nGrayLevels;
    idx2 = 0:nGrayLevels-1;
    for i = idx1
        pXplusY(i-1,k) = sum(tmp1(tmp1(:,1)==i,2));
    end
    
    for i = idx2 
        pXminusY(i+1,k) = sum(tmp2(tmp2(:,1)==i,2));
    end

    % These can be evaluated for all GLCMs simultaneously, no k-index
    % missing. We need the results further down so I keep it in the loop.
    out.sumAverage              = sum(bsxfun(@times,idx1',pXplusY));
    out.sumEntropy              = -nansum(pXplusY.*log(pXplusY)); %OK
    out.differenceEntropy       = -nansum(pXminusY.*log(pXminusY)); %OK
    out.differenceVariance(k)   = sum((idx2-out.dissimilarity(k)).^2'.*pXminusY(idx2+1,k)); %<----- N.B! Wrong implementation previously!! Dissimilarity is "difference Average"
    out.sumVariance(k)          = sum((idx1-out.sumAverage(k))'.^2.*pXplusY(idx1-1,k)); %<----- N.B! Wrong implementation previously AND in [1]
    
    HXY1(k)                     = -nansum(currentGLCM(sub)'.*log(pX(I,k).*pY(J,k))); %OK
    HXY2(k)                     = -nansum(pX(I,k).*pY(J,k).*log(pX(I,k).*pY(J,k))); %OK
    HX(k)                       = -nansum(pX(:,k).*log(pX(:,k))); %OK
    HY(k)                       = -nansum(pY(:,k).*log(pY(:,k))); %OK
    
    out.autoCorrelation(k)      = sum(I.*J.*currentGLCM(sub));
    out.clusterProminence(k)    = sum((I+J-uX(k)-uY(k)).^4.*currentGLCM(sub)); %OK
    out.clusterShade(k)         = sum((I+J-uX(k)-uY(k)).^3.*currentGLCM(sub)); %OK
    out.correlation(k)          = (out.autoCorrelation(k) - uX(k).*uY(k))./(sqrt(sX(k).*sY(k))); %OK
    
    out.informationMeasureOfCorrelation1(k) = (out.entropy(k)-HXY1(k))./(max(HX(k),HY(k))); %OK
    out.informationMeasureOfCorrelation2(k) = (1 - exp(-2.*(HXY2(k)-out.entropy(k))) ).^(1/2); %OK
    
end