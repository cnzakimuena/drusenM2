
function filterBW = logMed(otsuBW1, otsuBW2)

filterBW = zeros(size(otsuBW1));
for i = 1:size(otsuBW1,3)
    img1 = otsuBW1(:,:,i);
    img2 = otsuBW2(:,:,i);
    img1(img2) = 1; 
    img3 = medfilt2(img1);
    filterBW(:,:,i) = img3;
end

end