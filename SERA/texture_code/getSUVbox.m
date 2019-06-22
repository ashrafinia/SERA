function [SUVbox, maskBox] = getImgBox(volume,mask,isReSeg,ResegIntrval)

% COMPUTATION OF THE SMALLEST BOX CONTAINING THE ROI
[boxBound] = computeBoundingBox(mask);
maskBox = mask(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
SUVbox = volume(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

% Resegmentation
if nargin >2
    if isReSeg
        SUVbox(SUVbox<ResegIntrval(1)) = NaN;
        SUVbox(SUVbox>ResegIntrval(2)) = NaN;
    end
    
end
