function [GLroundedImg] = roundGL(Img , isGLrounding)
% This function rounds image intensity voxels to the nearest integer.

if isGLrounding
    GLroundedImg = round(Img);
else
    GLroundedImg = Img;
end 

end