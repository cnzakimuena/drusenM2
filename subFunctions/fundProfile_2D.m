
function areaProfile = fundProfile_2D(maxBW,rETDRS, sizeRed)

% RESULTS
%total drusen area
% 'radiusFac' is the conversion factor, 3000/1536 or 1.95 um/px 
radiusFac = 3000/1536*1/sizeRed; % um/px
totalArea = sum(maxBW(:))*radiusFac^2/1000^2;
aGR = zeros(1, size(rETDRS, 3));
for k = 1:size(rETDRS, 3)
    currArea = maxBW(logical(rETDRS(:,:,k)));
    realArea = sum(currArea(:))*radiusFac^2/1000^2;
    aGR(:, k) = realArea; % at given grid region, um^2
end
areaProfile = [totalArea aGR]; % um^2
% figure;imshow(maxBW,[])
% hold on
% plot(fovCenterX,fovCenterY,'.r')

% sub_maxBW = maxBW.*logical(regionsETDRS(:,:,1));
% figure;imshow(sub_maxBW,[])
% nnz(sub_maxBW) % sum of non-zeros for region 1 of chroidMap
% sum(sum((logical(regionsETDRS(:,:,1))))) % sum of non-zeros for ETDRS 1

end

