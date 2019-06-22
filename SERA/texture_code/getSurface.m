function surface = getSurface(mask,sizx,sizy,sizz)
area_vec = [sizy*sizz sizy*sizz...
    sizx*sizz sizx*sizz...
    sizx*sizy sizx*sizy]';%/1e2;%convert to cm^2
area_sum =0;
for idx1 = 2:size(mask,1)-1
    for idx2 = 2:size(mask,2)-1
        for idx3 = 2:size(mask,3)-1
            if mask(idx1,idx2,idx3)==1
                area_sum = area_sum + (~[mask(idx1-1,idx2,idx3),mask(idx1+1,idx2,idx3),mask(idx1,idx2-1,idx3),mask(idx1,idx2+1,idx3),mask(idx1,idx2,idx3-1),mask(idx1,idx2,idx3+1)])*area_vec;
            end
        end
    end
end
surface = area_sum;
