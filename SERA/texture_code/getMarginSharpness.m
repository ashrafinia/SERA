function MarginSharpness = getMarginSharpness(mask,image)
[Gx,Gy] = gradient(image);
G = sqrt(Gx.*Gx+Gy.*Gy);
grad_sum =0;
n = 0;
for idx1 = 2:size(mask,1)-1
    for idx2 = 2:size(mask,2)-1
        for idx3 = 2:size(mask,3)-1
            if mask(idx1,idx2,idx3)==1 && (mask(idx1-1,idx2,idx3)==0 || mask(idx1+1,idx2,idx3)==0 || mask(idx1,idx2-1,idx3)==0 || mask(idx1,idx2+1,idx3)==0 || mask(idx1,idx2,idx3-1)==0 || mask(idx1,idx2,idx3+1)==0);
               grad_sum = grad_sum + G(idx1,idx2,idx3);
               n = n+1;
            end
        end
    end
end
MarginSharpness = grad_sum./n;