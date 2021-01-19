
function jaccardBW = jaccardMap(filterBW, squareDim)

maxNum = size(filterBW, 3)-1;
jaccardBW0 = zeros(size(filterBW,1),size(filterBW,2), maxNum);
for i = 1:maxNum
    for ii = 1:size(jaccardBW0, 1)
        for iii = 1:size(jaccardBW0, 2)
            
            centerX = iii;
            centerY = ii;
 
            xmin = centerX-squareDim/2;
            ymin = centerY-squareDim/2;
            
            cropBW1 = logical(imcrop(filterBW(:,:,i),[xmin ymin ...
                squareDim squareDim]));
            cropBW2 = logical(imcrop(filterBW(:,:,i+1),[xmin ymin ...
                squareDim squareDim]));
            
            curr_jaccard = jaccard(cropBW1,cropBW2);
            if curr_jaccard > 0.5
                curr_jaccard = 1;
            else
                curr_jaccard = 0;
            end
            jaccardBW0(ii, iii, i) = curr_jaccard;
         
        end
    end
end

jaccardBW = sum(jaccardBW0, 3);
jaccardBW(jaccardBW <= 4) = 0; 
jaccardBW(jaccardBW > 4) = 1;
 
end