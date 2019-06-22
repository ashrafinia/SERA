function [metric] = getVolume_gt(ROIonlyPET,value,pixelW,sliceS)

metric=sum(ROIonlyPET(ROIonlyPET>=value))* pixelW * pixelW * sliceS;

end