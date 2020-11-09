%% mitometer
addpath(strcat(pwd,'/Simulation'))

extremaValues = [1,2,4,8,16,32];
noiseValues = [1,2,4,8,16,32,64];

numTrials = 3;
SNR = zeros(length(extremaValues),length(noiseValues),numTrials);
meanCompNum = zeros(length(extremaValues),length(noiseValues),numTrials);
medArea = zeros(length(extremaValues),length(noiseValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;

    noiseLevel = 0;
    for noiseTest = 1%:length(noiseValues)
        extremaDist = 0;
        noiseLevel = noiseLevel+1;
        for extremaTest = 1%:length(extremaValues)

            extremaDist = extremaDist+1;
            [mitoIm,mitoSim] = generateSpecificMitochondria(mitoLength,extremaValues(extremaTest),vert,255-noiseValues(noiseTest),100);

            comboMitoIm = zeros(100,100,5);

            tempIm1 = ones(100,100,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNR(extremaTest,noiseTest,trial) = signalVal/noiseVal;

            ImAll{1} = comboMitoIm;

            % User Parameters
            micronPerPixel = 0.138;
            minArea = 0.3;
            maxArea = 200;

            % User Image
            ImOG = ImAll{1};

            %remove diffuse bg and normalize image intensity
            [ImBgRemoved,ImMinMed,ImMedFiltered] = diffuseBgRemove(ImOG,(minArea)/micronPerPixel,(1)/micronPerPixel);

            %find optimal gaussian filter std Dev and hard threshold value
            [sig,thr,costs] = optimizeSigmaThresh(ImBgRemoved);

            %gaussian filter for low-pass filter and smooth discontinuities.
            ImGaussFiltered = gaussFilter(ImBgRemoved,sig);

            %threshold image to remove adjacent mitochondria connections and bg.
            ImMask = thresholdImage(ImGaussFiltered,thr);

            %apply mask
            ImMaskThresholded = areaThreshold(ImMask,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));

            for k = 1:size(ImMaskThresholded,3)
                ImMaskThresholded(:,:,k) = imclearborder(ImMaskThresholded(:,:,k));
            end

            Im = ImMaskThresholded.*ImOG;

            numComps = zeros(1,size(ImMaskThresholded,3));
            for i = 1:size(ImMaskThresholded,3)
                compStart = bwconncomp(ImMaskThresholded(:,:,i));
                numComps(i) = compStart.NumObjects;
            end
                        
            areaList = zeros(1,length(compStart.PixelIdxList));
            for i = 1:length(compStart.PixelIdxList)
                areaList(i) = length(compStart.PixelIdxList{i});
            end
            
            medAreaTemp = median(areaList);
            medArea(extremaTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNum(extremaTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,extremaTest;meanCompNumTemp,noiseValues(noiseTest),extremaValues(extremaTest)]
            meanCompNum(:,:,trial)
            trial
            1

        end
    end
end
%%
% small med filt + gauss,thr algorithm

extremaValues = [1,2,4,8,16,32];
noiseValues = [1,2,4,8,16,32,64];

numTrials = 3;
SNRMF3 = zeros(length(extremaValues),length(noiseValues),numTrials);
meanCompNumMF3 = zeros(length(extremaValues),length(noiseValues),numTrials);
medAreaMF3 = zeros(length(extremaValues),length(noiseValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;

    noiseLevel = 0;
    for noiseTest = 1:length(noiseValues)
        extremaDist = 0;
        noiseLevel = noiseLevel+1;
        for extremaTest = 1:length(extremaValues)

            extremaDist = extremaDist+1;
            [mitoIm,mitoSim] = generateSpecificMitochondria(mitoLength,extremaValues(extremaTest),vert,255-noiseValues(noiseTest),100);

            comboMitoIm = zeros(100,100,5);

            tempIm1 = ones(100,100,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNRMF3(extremaTest,noiseTest,trial) = signalVal/noiseVal;

            ImAll{1} = comboMitoIm;

            % User Parameters
            micronPerPixel = 0.138;
            minArea = 0.3;
            maxArea = 200;

            % User Image
            ImOG = ImAll{1};

            %one-sized medfilt 3x3
            ImBgRemoved = zeros(size(comboMitoIm));
            for frameNum = 1:size(comboMitoIm,3)
                ImBgRemoved(:,:,frameNum) = comboMitoIm(:,:,frameNum)-medfilt2(comboMitoIm(:,:,frameNum),[3,3]);
            end
            
            %find optimal gaussian filter std Dev and hard threshold value
            [sig,thr,costs] = optimizeSigmaThresh(ImBgRemoved);

            %gaussian filter for low-pass filter and smooth discontinuities.
            ImGaussFiltered = gaussFilter(ImBgRemoved,sig);

            %threshold image to remove adjacent mitochondria connections and bg.
            ImMask = thresholdImage(ImGaussFiltered,thr);

            %apply mask
            ImMaskThresholded = areaThreshold(ImMask,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));

            for k = 1:size(ImMaskThresholded,3)
                ImMaskThresholded(:,:,k) = imclearborder(ImMaskThresholded(:,:,k));
            end

            Im = ImMaskThresholded.*ImOG;

            numComps = zeros(1,size(ImMaskThresholded,3));
            for i = 1:size(ImMaskThresholded,3)
                compStart = bwconncomp(ImMaskThresholded(:,:,i));
                numComps(i) = compStart.NumObjects;
            end
                        
            areaList = zeros(1,length(compStart.PixelIdxList));
            for i = 1:length(compStart.PixelIdxList)
                areaList(i) = length(compStart.PixelIdxList{i});
            end
            
            medAreaTemp = median(areaList);
            medAreaMF3(extremaTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNumMF3(extremaTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,extremaTest;meanCompNumTemp,noiseValues(noiseTest),extremaValues(extremaTest)]
            meanCompNumMF3(:,:,trial)
            trial
            2

        end
    end
end

% big med filt + gauss,thr algorithm

extremaValues = [1,2,4,8,16,32];
noiseValues = [1,2,4,8,16,32,64];

numTrials = 3;
SNRMF14 = zeros(length(extremaValues),length(noiseValues),numTrials);
meanCompNumMF14 = zeros(length(extremaValues),length(noiseValues),numTrials);
medAreaMF14 = zeros(length(extremaValues),length(noiseValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;

    noiseLevel = 0;
    for noiseTest = 1:length(noiseValues)
        extremaDist = 0;
        noiseLevel = noiseLevel+1;
        for extremaTest = 1:length(extremaValues)

            extremaDist = extremaDist+1;
            [mitoIm,mitoSim] = generateSpecificMitochondria(mitoLength,extremaValues(extremaTest),vert,255-noiseValues(noiseTest),100);

            comboMitoIm = zeros(100,100,5);

            tempIm1 = ones(100,100,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNRMF14(extremaTest,noiseTest,trial) = signalVal/noiseVal;

            ImAll{1} = comboMitoIm;

            % User Parameters
            micronPerPixel = 0.138;
            minArea = 0.3;
            maxArea = 200;

            % User Image
            ImOG = ImAll{1};

            %one-sized medfilt 3x3
            ImBgRemoved = zeros(size(comboMitoIm));
            for frameNum = 1:size(comboMitoIm,3)
                ImBgRemoved(:,:,frameNum) = comboMitoIm(:,:,frameNum)-medfilt2(comboMitoIm(:,:,frameNum),[14,14]);
            end
            
            %find optimal gaussian filter std Dev and hard threshold value
            [sig,thr,costs] = optimizeSigmaThresh(ImBgRemoved);

            %gaussian filter for low-pass filter and smooth discontinuities.
            ImGaussFiltered = gaussFilter(ImBgRemoved,sig);

            %threshold image to remove adjacent mitochondria connections and bg.
            ImMask = thresholdImage(ImGaussFiltered,thr);

            %apply mask
            ImMaskThresholded = areaThreshold(ImMask,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));

            for k = 1:size(ImMaskThresholded,3)
                ImMaskThresholded(:,:,k) = imclearborder(ImMaskThresholded(:,:,k));
            end

            Im = ImMaskThresholded.*ImOG;

            numComps = zeros(1,size(ImMaskThresholded,3));
            for i = 1:size(ImMaskThresholded,3)
                compStart = bwconncomp(ImMaskThresholded(:,:,i));
                numComps(i) = compStart.NumObjects;
            end
                        
            areaList = zeros(1,length(compStart.PixelIdxList));
            for i = 1:length(compStart.PixelIdxList)
                areaList(i) = length(compStart.PixelIdxList{i});
            end
            
            medAreaTemp = median(areaList);
            medAreaMF14(extremaTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNumMF14(extremaTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,extremaTest;meanCompNumTemp,noiseValues(noiseTest),extremaValues(extremaTest)]
            meanCompNumMF14(:,:,trial)
            trial
            3

        end
    end
end

% medfilt algorithm + small gauss + thr algorithm


extremaValues = [1,2,4,8,16,32];
noiseValues = [1,2,4,8,16,32,64];

numTrials = 3;
SNRGF3 = zeros(length(extremaValues),length(noiseValues),numTrials);
meanCompNumGF3 = zeros(length(extremaValues),length(noiseValues),numTrials);
medAreaGF3 = zeros(length(extremaValues),length(noiseValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;

    noiseLevel = 0;
    for noiseTest = 1:length(noiseValues)
        extremaDist = 0;
        noiseLevel = noiseLevel+1;
        for extremaTest = 1:length(extremaValues)

            extremaDist = extremaDist+1;
            [mitoIm,mitoSim] = generateSpecificMitochondria(mitoLength,extremaValues(extremaTest),vert,255-noiseValues(noiseTest),100);

            comboMitoIm = zeros(100,100,5);

            tempIm1 = ones(100,100,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNRGF3(extremaTest,noiseTest,trial) = signalVal/noiseVal;

            ImAll{1} = comboMitoIm;

            % User Parameters
            micronPerPixel = 0.138;
            minArea = 0.3;
            maxArea = 200;

            % User Image
            ImOG = ImAll{1};
            
            %remove diffuse bg and normalize image intensity
            [ImBgRemoved,ImMinMed,ImMedFiltered] = diffuseBgRemove(ImOG,(minArea)/micronPerPixel,(1)/micronPerPixel);

            %find optimal gaussian filter std Dev and hard threshold value
            [sig,thr,costs] = optimizeOnlyThresh(ImBgRemoved,0.33);

            %gaussian filter for low-pass filter and smooth discontinuities.
            ImGaussFiltered = gaussFilter(ImBgRemoved,sig);

            %threshold image to remove adjacent mitochondria connections and bg.
            ImMask = thresholdImage(ImGaussFiltered,thr);

            %apply mask
            ImMaskThresholded = areaThreshold(ImMask,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));

            for k = 1:size(ImMaskThresholded,3)
                ImMaskThresholded(:,:,k) = imclearborder(ImMaskThresholded(:,:,k));
            end

            Im = ImMaskThresholded.*ImOG;

            numComps = zeros(1,size(ImMaskThresholded,3));
            for i = 1:size(ImMaskThresholded,3)
                compStart = bwconncomp(ImMaskThresholded(:,:,i));
                numComps(i) = compStart.NumObjects;
            end
                        
            areaList = zeros(1,length(compStart.PixelIdxList));
            for i = 1:length(compStart.PixelIdxList)
                areaList(i) = length(compStart.PixelIdxList{i});
            end
            
            medAreaTemp = median(areaList);
            medAreaGF3(extremaTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNumGF3(extremaTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,extremaTest;meanCompNumTemp,noiseValues(noiseTest),extremaValues(extremaTest)]
            meanCompNumGF3(:,:,trial)
            trial
            4

        end
    end
end

% medfilt algorithm + small gauss + thr algorithm


extremaValues = [1,2,4,8,16,32];
noiseValues = [1,2,4,8,16,32,64];

numTrials = 3;
SNRGF5 = zeros(length(extremaValues),length(noiseValues),numTrials);
meanCompNumGF5 = zeros(length(extremaValues),length(noiseValues),numTrials);
medAreaGF5 = zeros(length(extremaValues),length(noiseValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;

    noiseLevel = 0;
    for noiseTest = 1:length(noiseValues)
        extremaDist = 0;
        noiseLevel = noiseLevel+1;
        for extremaTest = 1:length(extremaValues)

            extremaDist = extremaDist+1;
            [mitoIm,mitoSim] = generateSpecificMitochondria(mitoLength,extremaValues(extremaTest),vert,255-noiseValues(noiseTest),100);

            comboMitoIm = zeros(100,100,5);

            tempIm1 = ones(100,100,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNRGF5(extremaTest,noiseTest,trial) = signalVal/noiseVal;

            ImAll{1} = comboMitoIm;

            % User Parameters
            micronPerPixel = 0.138;
            minArea = 0.3;
            maxArea = 200;

            % User Image
            ImOG = ImAll{1};

            %remove diffuse bg and normalize image intensity
            [ImBgRemoved,ImMinMed,ImMedFiltered] = diffuseBgRemove(ImOG,(minArea)/micronPerPixel,(1)/micronPerPixel);

            %find optimal gaussian filter std Dev and hard threshold value
            [sig,thr,costs] = optimizeOnlyThresh(ImBgRemoved,0.5);

            %gaussian filter for low-pass filter and smooth discontinuities.
            ImGaussFiltered = gaussFilter(ImBgRemoved,sig);

            %threshold image to remove adjacent mitochondria connections and bg.
            ImMask = thresholdImage(ImGaussFiltered,thr);

            %apply mask
            ImMaskThresholded = areaThreshold(ImMask,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));

            for k = 1:size(ImMaskThresholded,3)
                ImMaskThresholded(:,:,k) = imclearborder(ImMaskThresholded(:,:,k));
            end

            Im = ImMaskThresholded.*ImOG;

            numComps = zeros(1,size(ImMaskThresholded,3));
            for i = 1:size(ImMaskThresholded,3)
                compStart = bwconncomp(ImMaskThresholded(:,:,i));
                numComps(i) = compStart.NumObjects;
            end
                        
            areaList = zeros(1,length(compStart.PixelIdxList));
            for i = 1:length(compStart.PixelIdxList)
                areaList(i) = length(compStart.PixelIdxList{i});
            end
            
            medAreaTemp = median(areaList);
            medAreaGF5(extremaTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNumGF5(extremaTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,extremaTest;meanCompNumTemp,noiseValues(noiseTest),extremaValues(extremaTest)]
            meanCompNumGF5(:,:,trial)
            trial
            5
        end
    end
end