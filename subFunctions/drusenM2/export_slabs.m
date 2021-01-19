function export_slabs(inner_slab, outer_slab, destFolder, folderName)

for i = 1:size(inner_slab, 3)
    A = inner_slab(:,:,i);
    imwrite(A,fullfile([destFolder,'\Results\', folderName,'\innerImage_' num2str(i,'%02.f') '.png']));
%     B = inner_slab(:,:,size(inner_slab, 3) - (i-1));
%     imwrite(B,[ 'innerImage_' num2str(size(inner_slab, 3)+i,'%02.f') '.png']);
end
for ii = 1:size(outer_slab, 3)
    A2 = outer_slab(:,:,ii);
    imwrite(A2,fullfile([destFolder,'\Results\', folderName, '\outerImage_' num2str(ii,'%02.f') '.png']));
%     B2 = outer_slab(:,:,size(outer_slab, 3) - (ii-1));
%     imwrite(B2,[ 'outerImage_' num2str(size(outer_slab, 3)+ii,'%02.f') '.png']);    
end

end