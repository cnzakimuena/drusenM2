
function call_drusenM2()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name - drusenM2
% Creation Date - 17th January 2021
% Author - Charles Belanger Nzakimuena
% Website - https://www.ibis-space.com/
%
% Description - 
%   The 'drusenM2' algorithm uses Bruch’s membrane (BM) segemntation data 
%   towards drusen detection. Area values corresponding to ETDRS subfieds 
%   are provided in table format.
%
% Example -
%		call_drusenM2()
%
% License - MIT
%
% Change History -
%                   17th January 2021 - Creation by Charles Belanger Nzakimuena
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('./subfunctions'))

%% list names of folders inside the patients folder

currentFolder = pwd;
patientsFolder = fullfile(currentFolder, 'processed');
myDir = dir(patientsFolder);
dirFlags = [myDir.isdir] & ~strcmp({myDir.name},'.') & ~strcmp({myDir.name},'..');
nameFolds = myDir(dirFlags);

%% for each 3x3 subfolder, turn segmented data into network graph

% get table row count
rowCount = 0;
for g = 1:numel(nameFolds)
    folder2 = fullfile(patientsFolder, nameFolds(g).name);
    patientDir2 = dir(fullfile(folder2, 'Results'));
    dirFlags2 = [patientDir2.isdir] & ~strcmp({patientDir2.name},'.') & ~strcmp({patientDir2.name},'..');
    subFolders2 = patientDir2(dirFlags2);
    rowCount = rowCount + numel(subFolders2);
end

col = zeros(rowCount,1);
colc = cell(rowCount,1);
% aTable2 = table(colc,col,col,col,col,col,col,col,col,col,col,col,...
%     'VariableNames',{'id' 'totalArea' 'region1' 'region2' 'region3' 'region4'...
%     'region5' 'contrastRatio' 'jaccard' 'dice' 'sensitivity' 'specificity'});
aTable2 = table(colc,col,col,col,col,col,col,col,...
    'VariableNames',{'id' 'totalArea' 'region1' 'region2' 'region3' 'region4'...
    'region5' 'contrastRatio'});

tableRow = 0;
for i = 1:numel(nameFolds)
    
    % assemble patient folder string
    folder = fullfile(patientsFolder, nameFolds(i).name);
    
    % add line to LOG
    disp(logit(folder, ['Initiating drusenM2; ' nameFolds(i).name ' folder']))
    
    patientDir = dir(fullfile(folder, 'Results'));
    dirFlags = [patientDir.isdir] & ~strcmp({patientDir.name},'.') & ~strcmp({patientDir.name},'..');
    subFolders = patientDir(dirFlags);

    for k = 1:numel(subFolders)
        
        nameFold = subFolders(k).name;
        scanType = nameFold(1:2);
        
        if strcmp(scanType, '3m')
            
            load(fullfile(folder,'Results', nameFold,'segmentation.mat'));
            load(fullfile(folder,'Results', nameFold,'scanInfo.mat'));
            load(fullfile(folder,'Results', nameFold, 'ETDRS_grid','2DregionsETDRS.mat'));
            load(fullfile(folder,'Results', nameFold, 'ETDRS_grid','3DregionsETDRS.mat'));
            sizeRed = scanTag{2};
            
            clear RPEb RPEt RVIf volumeFlow
            
            %% PRE-PROCESSING STEPS 3, 4 : INNER AND OUTER SLABS
            disp('begin STEPS 3, 4 : inner and outer slabs')
            % Two slabs each containing 14 C-scans were extracted.
            % The inner and outer slab are located at 42-56 pixels and 0-14 pixels
            % above the Bruch’s membrane respectively.
            % Slab bounderies were located through trial and error
            % INNER SLAB ROI EXTRACTION
            % slab sampled as follows
            innerLoc_t = 56; % px above BM (1.953um/px)
            innerLoc_b = 42; % px above BM (1.953um/px)
            [innerSlab, innerSurf_t, innerSurf_b] = volROI2(volumeStruc, lBM, innerLoc_t, innerLoc_b);
            clear innerLoc_t innerLoc_b
            % figure;imshow3D(innerSlab,[])
            % figure;imshow3D(volumeStruc,[],'plot',cat(3,lBM, innerSurf_t, innerSurf_b),'LineWidth',2)
            % OUTER SLAB ROI EXTRACTION
            % slab sampled as follows
            outerLoc_t = 18; % px above BM (1.953um/px)
            outerLoc_b = 4; % px above BM (1.953um/px)
            [outerSlab, outerSurf_t, outerSurf_b] = volROI2(volumeStruc, lBM, outerLoc_t, outerLoc_b);
            clear outerLoc_t outerLoc_b
%             clear volumeStruc lBM
            % figure;imshow3D(outerSlab,[])
            % figure;imshow3D(volumeStruc,[],'plot',cat(3,lBM, outerSurf_t, outerSurf_b),'LineWidth',2)
            
%             chartColors1.c2 = rgb('FireBrick');
%             chartColors2.c2 = rgb('LightPink');
%             figure;imshow(volumeStruc(:,:,155),[])
%             hold on;
%             plot(innerSurf_t(155,:),'--','color',chartColors2.c2,'LineWidth',1.5);
%             plot(innerSurf_b(155,:),'--','color',chartColors2.c2,'LineWidth',1.5);
%             plot(outerSurf_t(155,:),'--','color',chartColors1.c2,'LineWidth',1.5);
%             plot(outerSurf_b(155,:),'--','color',chartColors1.c2,'LineWidth',1.5);
             
            % % inner and outer slabs images export
            export_slabs(innerSlab, outerSlab, folder, nameFold)

            disp('end STEPS 3, 4 : inner and outer slabs')
            
            %% PRE-PROCESSING STEP 5 : HYBRID CONTRAST RATIO
            disp('begin STEPS 5 : hybrid contrast ratio')
            % slabs c-scan reference study size reduction
            scaleAdjust = 300/1536;
            
            innerSlab_sm = imresize3(innerSlab,[round(size(innerSlab,1)*scaleAdjust) ...
                round(size(innerSlab,2)*scaleAdjust) size(innerSlab,3)]);
            outerSlab_sm = imresize3(outerSlab,[round(size(outerSlab,1)*scaleAdjust) ...
                round(size(outerSlab,2)*scaleAdjust) size(outerSlab,3)]);
            clear scaleAdjust innerSlab innerSurf_t innerSurf_b
            clear outerSlab outerSurf_t outerSurf_b
            % inner GCxy, LCxy and hybrid contast map calculation
            GCxy_inner = GCxy_calc(innerSlab_sm);
            LCxy_inner = LCxy_calc(innerSlab_sm);
            hybridCM_inner = LCxy_inner.*GCxy_inner;
            % figure;imshow([GCxy_inner LCxy_inner hybridCM_inner],[])
            % outer GCxy, LCxy and hybrid contast map calculation
            GCxy_outer = GCxy_calc(outerSlab_sm);
            LCxy_outer = LCxy_calc(outerSlab_sm);
            hybridCM_outer = LCxy_outer.*GCxy_outer;
            % figure;imshow([GCxy_outer LCxy_outer hybridCM_outer],[])
            contrast_ratio = mean([max(max(hybridCM_inner)) max(max(hybridCM_outer))]);
            clear GCxy_inner LCxy_inner hybridCM_inner
            clear GCxy_outer LCxy_outer hybridCM_outer
            disp('end STEPS 5 : hybrid contrast ratio')
            
            %% DRUSEN SEGMENTATION STEP 6, 7 : AVERAGED C-SCANS, INVERSED OUTER SLAB
            disp('begin STEPS 6, 7 : averaged C-scans, inversed outer slab')
            dodecuplet_inner = cscan3_mean(innerSlab_sm);
            dodecuplet_outer = cscan3_mean(outerSlab_sm);
            inversed_outer = imcomplement(dodecuplet_outer);
            clear innerSlab_sm outerSlab_sm
            disp('end STEPS 6, 7 : averaged C-scans, inversed outer slab')
            
            %% DRUSEN SEGMENTATION STEP 8, 9 : BINARY IMAGES
            disp('begin STEPS 8, 9 : binary images')
            % inner slab dual Top-Hat transform and Otsu implementation
            otsuBW_inner1 = logical(tophat_Otsu1(dodecuplet_inner, 60));
            otsuBW_inner2 = logical(tophat_Otsu2(dodecuplet_inner, 15));
%             clear dodecuplet_inner

            % figure; imshow3D(otsuBW_inner1, [])
            % figure; imshow3D(otsuBW_inner2, [])
            % outer slab dual Top-Hat transform and Otsu implementation
            otsuBW_outer1 = logical(tophat_Otsu3(inversed_outer, 60));
            otsuBW_outer2 = logical(tophat_Otsu4(inversed_outer, 15));
%             clear dodecuplet_outer inversed_outer
            % figure; imshow3D(otsuBW_outer1, [])
            % figure; imshow3D(otsuBW_outer2, [])

            % *ADDED MORPHOLOGICAL OPERATIONS START*
            img1 = zeros(size(otsuBW_inner1));
            img2 = zeros(size(otsuBW_inner1));
            img3 = zeros(size(otsuBW_inner1));
            img4 = zeros(size(otsuBW_inner1));
            for kk = 1:size(otsuBW_inner1,3)
                
                currImg3 = bwareaopen(otsuBW_inner1(:,:,kk), 120); % default 20
                if kk <= 7
                    se1 = strel('disk',1);
                    currImg2 = imerode(otsuBW_inner1(:,:,kk),se1);
                    currImg3 = bwareaopen(currImg2, 120); % default 20
                end
                img1(:,:,kk) = currImg3;
                
                currImg5 = bwareaopen(otsuBW_inner2(:,:,kk), 180); % default 40
                if kk <= 7
                    se2 = strel('disk',1);
                    currImg4 = imerode(otsuBW_inner2(:,:,kk),se2);
                    currImg5 = bwareaopen(currImg4, 180); % default 20
                end
                img2(:,:,kk) = currImg5;

                se3 = strel('disk',2);
                currImg6 = imerode(otsuBW_outer1(:,:,kk),se3);
                currImg7 = bwareaopen(currImg6, 240); % default 40
                img3(:,:,kk) = currImg7;
                
%                 se4 = strel('disk',2);
%                 currImg8 = imerode(otsuBW_outer2(:,:,i),se4);
                currImg9 = bwareaopen(otsuBW_outer2(:,:,kk), 180); % default 20
                img4(:,:,kk) = currImg9;
                
                imgTop = [dodecuplet_inner(:,:,kk) dodecuplet_inner(:,:,kk) inversed_outer(:,:,kk) inversed_outer(:,:,kk)];
                imgBottom = [currImg3 currImg5 currImg7 currImg9];
                imgT(:,:,kk) = [imgTop; imgBottom];
                
            end
            img1 = logical(img1);
            img2 = logical(img2);
            img3 = logical(img3);
            img4 = logical(img4);
%             figure; imshow3D(imgT,[])
            % *ADDED MORPHOLOGICAL OPERATIONS END*

            % inner slab logicalOR and median filtering
            filterBW_inner = logMed(img1, img2);
%             clear otsuBW_inner1 otsuBW_inner2
            % figure; imshow3D(filterBW_inner, [])
            % outer slab logicalOR and median filtering
            filterBW_outer =logMed(img3, img4);
%             clear otsuBW_outer1 otsuBW_outer2
            % figure; imshow3D(filterBW_outer, [])
            disp('end STEPS 8, 9 : binary images')
            
            %% POST-PROCESSING STEP 10, 11, 12 : CONTINUITY, OPENING, LOGICAL OR
            disp('begin STEPS 10, 11, 12 : continuity, opening, logical OR')
            square_Dim = 3; % real desired window dimension in px
            jaccard_inner = logical(jaccardMap(filterBW_inner, square_Dim));
            clear filterBW_inner
            % figure;imshow(jaccard_inner,[])
            jaccard_outer = logical(jaccardMap(filterBW_outer, square_Dim));
            clear filterBW_outer square_Dim    
            
            % figure;imshow(jaccard_outer2,[])
            drusenMask = logMed(jaccard_inner, jaccard_outer);
            clear jaccard_inner jaccard_outer2
            % figure;imshow(drusenMask,[])
            disp('end STEPS 10, 11, 12 : continuity, opening, logical OR')

            numrows2 = round(size(drusenMask, 1)*1536/300*sizeRed);
            numcols2 = round(size(drusenMask, 2)*1536/300*sizeRed);
            drusenMask_BW = imbinarize(imresize(drusenMask,[numrows2 numcols2]));
            imwrite(drusenMask_BW,fullfile([folder,'\Results\', nameFold, '\drusenMaskM2_BW' '.png']));
            %figure;imshow(drusenMask_BW,[])               

%             % obtain validation parameters
%             if strcmp(patientsType,'AMD')
%                 valid_folder = 'C:\Users\...';
%                 valid_Mask = imread(fullfile([valid_folder,'\',nameFolds(i).name,'\',nameFold, '.png']));
%                 if size(valid_Mask,3)==3
%                     valid_Mask = rgb2gray(valid_Mask);
%                 end
%                 level = graythresh(valid_Mask);
%                 valid_MaskBW = imbinarize(valid_Mask,level);
%                 valid_MaskBW = imresize(valid_MaskBW,[numrows2 numcols2]);
%             elseif strcmp(patientsType,'normal') || strcmp(patientsType,'normal_U50') 
%                 valid_MaskBW = imbinarize(zeros(size(drusenMask_BW)));
%             end
%             jaccard_coef = jaccard(valid_MaskBW,drusenMask_BW);
%             dice_coef = dice(valid_MaskBW,drusenMask_BW);
%             TP=0;FP=0;TN=0;FN=0;
%             for w=1:size(drusenMask_BW, 1)
%                 for ww=1:400
%                     if(valid_MaskBW(w,ww)==1 && drusenMask_BW(w,ww)==1)
%                         TP=TP+1;
%                     elseif(valid_MaskBW(w,ww)==0 && drusenMask_BW(w,ww)==1)
%                         FP=FP+1;
%                     elseif(valid_MaskBW(w,ww)==0 && drusenMask_BW(w,ww)==0)
%                         TN=TN+1;
%                     else
%                         FN=FN+1;
%                     end
%                 end
%             end
%             sensitivity = TP/(TP+FN);
%             specificity = TN/(TN+FP);

            % 2D ETDRS regions setup
            disp('begin fundProfile_2D')
            aProfile = fundProfile_2D(drusenMask_BW, regionsETDRS, sizeRed);
            disp('end fundProfile_2D')
            
%             aProfile = [aProfile contrast_ratio jaccard_coef dice_coef sensitivity specificity];
            aProfile = [aProfile contrast_ratio];
            
            % For left eye, ETDRS regions must be modified from OD nomenclature
            % to OS nomenclature
            if contains(nameFold, '_OS_')
                aRegion3 = aProfile(6);
                aRegion5 = aProfile(4);
                aProfile(4) = aRegion3;
                aProfile(6) = aRegion5;
                
            end
            
            tableRow = tableRow + 1;
            
            aTable2{tableRow,'id'} = {nameFold};
            aTable2{tableRow,'totalArea'} = aProfile(1);
            aTable2{tableRow,'region1'}  = aProfile(2);
            aTable2{tableRow,'region2'} = aProfile(3);
            aTable2{tableRow,'region3'} = aProfile(4);
            aTable2{tableRow,'region4'} = aProfile(5);
            aTable2{tableRow,'region5'} = aProfile(6);
            aTable2{tableRow,'contrastRatio'} = aProfile(7);
%             aTable2{tableRow,'jaccard'} = aProfile(8);
%             aTable2{tableRow,'dice'} = aProfile(9);
%             aTable2{tableRow,'sensitivity'} = aProfile(10);
%             aTable2{tableRow,'specificity'} = aProfile(11);
            
        end
    end
    
end

fileName1 = fullfile(patientsFolder,'aTable2.xls');
writetable(aTable2,fileName1)

disp(logit(patientsFolder,'Done drusenM2'))

            
            