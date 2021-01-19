function LCxy = LCxy_calc(slab)

avgIm = zeros(size(slab(:,:,1)));
for i = 1:size(avgIm, 1)
    for ii = 1:size(avgIm, 2)
        avgIm(i,ii) = mean(slab(i, ii, :));
    end
end
normIm = (avgIm-min(min(avgIm)))/(max(max(avgIm))-min(min(avgIm)));

res_cscan = 3000/round(size(slab,1)); % fundus resolution um/px
wc_diaDim = 60; % real desired diameter dimension in um
ws_diaDim = 130; % real desired diameter dimension in um
radiusFac = 1/res_cscan; % conversion factor px/um
wc_radius = round((wc_diaDim*radiusFac)/2);
ws_radius = round((ws_diaDim*radiusFac)/2);
cscanDimX = size(normIm,2); % actual c-scan X dimension in um
cscanDimY = size(normIm,1); % actual c-scan Y dimension in um
cscanSizeX = round(cscanDimX);
cscanSizeY = round(cscanDimY);
[columnsInCscan, rowsInCscan] = meshgrid(1:cscanSizeX, 1:cscanSizeY);
wc = zeros(size(normIm));
ws = zeros(size(normIm));
for i = 1:size(normIm, 1)
    for ii = 1:size(normIm, 2)
        centerX = ii;
        centerY = i;
        wc_circlePixels = (rowsInCscan - centerY).^2 + (columnsInCscan - centerX).^2 <= wc_radius.^2;
        % figure;imshow(wc_circlePixels,[])
        ws_circlePixels = (rowsInCscan - centerY).^2 + (columnsInCscan - centerX).^2 <= ws_radius.^2;
        ws_circlePixels(wc_circlePixels) = 0;
        % figure;imshow(ws_circlePixels,[])
        wc(i, ii) = mean(normIm(logical(wc_circlePixels)));
        ws(i, ii) = mean(normIm(logical(ws_circlePixels)));
    end
end
LCxy = wc./ws;

end
