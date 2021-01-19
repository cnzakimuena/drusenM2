function GCxy = GCxy_calc(slab)

avgIm = zeros(size(slab(:,:,1)));
for i = 1:size(avgIm, 1)
    for ii = 1:size(avgIm, 2)
        avgIm(i,ii) = mean(slab(i, ii, :));
    end
end
normIm = (avgIm-min(min(avgIm)))/(max(max(avgIm))-min(min(avgIm)));
GCxy = normIm/mean(mean(normIm));

end