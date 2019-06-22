function [outROI] = getMutualROI(ROI1, ROI2)
% This function receives two ROIs as a result of two resegmentation methods
% (range reseg. and outlier reseg.) and return a third ROI containing the
% mutual voxels inside both ROIs. 


tmp1 = ROI1.* ROI2; 
tmp2=tmp1; 
tmp2(~isnan(tmp1)) = 1; 
outROI = tmp2 .* ROI1;