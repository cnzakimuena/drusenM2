function meanGroups = cscan3_mean(slab)

groupSize = size(slab, 3)-2;
meanGroups = zeros(size(slab,1),size(slab,2), groupSize);
for i = 1:groupSize
       sum_cscans = slab(:,:,i)+slab(:,:,i+1)+slab(:,:,i+2);
       avg_3cscans = sum_cscans/3;
       normIm = (avg_3cscans-min(min(avg_3cscans)))/(max(max(avg_3cscans))-min(min(avg_3cscans)));
       meanGroups(:,:,i) = normIm;
end

end