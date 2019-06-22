function EIG = getPCAeig(SUVbox, pixelW,sliceS)

% This function calculates eigenvalues of PCA for an ROI


[x,y,z]=ind2sub(size(SUVbox),find(~isnan(SUVbox)));
[~,EIG]=pcafun([x*pixelW, y*pixelW, z*sliceS]); 

end