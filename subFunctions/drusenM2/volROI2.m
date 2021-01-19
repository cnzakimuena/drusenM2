
function [slab_enface, slab_t, slab_b] = volROI2(Volume1, l_lim, wt, wb)

% INPUT DESCRIPTION
% 'Volume1'
%   The scan volume, with A-scan depth of 1536 pixels, B-scan width of
%   1527 pixels, and 300 B-scans in the volume direction
% 'l_lim'
%   The # rows corresponds to # of B-scans, and the value of each column
%   in a given row constitutes each Bruch membrane segmentation line
% figure;imshow3D(Volume1,[],'plot',l_lim,'LineWidth',2);

% OUPUT DESCRIPTION
% 'innerSlab'
% 'inner_t', 'inner_b'

slab_t = l_lim - wt;
slab_b = l_lim - wb;
% figure;imshow3D(vol1,[],'plot',cat(3,l_lim, inner_t, inner_b),'LineWidth',2)
% figure;imshow3D(cropVol1,[],'plot',cat(3,l_lim, inner_t, inner_b),'LineWidth',2)
% figure;imshow(Volume1(:,:,150),[])
% hold on;
% plot(inner_t(150,:),'b--','LineWidth',1);
% plot(inner_b(150,:),'b--','LineWidth',1);

slabThickness = abs(wt - wb);
innerSlab = zeros([abs(wt - wb) size(l_lim, 2) size(l_lim, 1)]);
for i = 1:slabThickness
    % begins with C-scans furthest away from BM
    currC_scan = l_lim - (wt - (i - 1));
    for k = 1:size(currC_scan, 1)
        for kk = 1:size(currC_scan, 2)
            innerSlab(i, kk, k) = Volume1(currC_scan(k,kk), kk, k);
        end
    end
end
curr_enface = zeros([size(innerSlab,3) size(innerSlab,2)]);
slab_enface = zeros([size(innerSlab,3) size(innerSlab,2) size(innerSlab,1)]);
for ii = 1:size(innerSlab,1)
    for g = 1:size(innerSlab,3)
        curr_row = innerSlab(ii,:,g);
        curr_enface(g,:) = flip(curr_row);
    end
    slab_enface(:,:,ii) = curr_enface;
end
slab_enface = imresize3(slab_enface,[1536 size(Volume1, 2) slabThickness]);
% figure; imshow3D(slab_enface,[])

% *uncomment below to obtain cropped ROI*
% vol1 = zeros(size(Volume1));
% cropVol1 = zeros(size(Volume1));
% for vv = 1:size(Volume1, 3) % for all B-scans
%     im1 = Volume1(:,:,vv); % current B-scan
%     
%     % create a mask with ones in the area of interest
%     cropMask = zeros(size(im1));
%     
%     % for the current B-scan, iterate through each A-scan column
%     for rr = 1:size(cropMask, 2)
%         % assign 1 at given A-scan slab location
%         cropMask(l_lim(vv,rr)-wt:l_lim(vv,rr)-wb,rr) = 1;
%     end
%    
%     im1 = cropMask.*im1;
%     vol1(:,:,vv) = im1; % assign cropped image to new volume
%     cropVol1(:,:,vv) = cropMask;
% end

% *uncomment below to obtain volume in the enface direction*
% vol1_enface = imrotate3(imrotate3(Volume1, -90, [1 0 0]), 180, [0 0 1]);
% figure; imshow3D(vol1_enface, []);

end
