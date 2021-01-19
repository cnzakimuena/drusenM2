
function img1 = tophat_Otsu1(img0, seNum)

% structuring element with disk of given radius in px
se = strel('disk',seNum); 
% slab Top-Hat transform implementation
tophatFiltered = imtophat(img0,se);
% figure; imshow3D(tophatFiltered, [])

% slab Otsu thresholding
img1 = zeros(size(tophatFiltered));
for i = 1:size(tophatFiltered,3)
%     level = graythresh(tophatFiltered(:,:,i));
    level = 0.3250;
    img1(:,:,i) = imbinarize(tophatFiltered(:,:,i), level);
end
% figure; imshow3D(img1, [])

end