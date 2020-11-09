%% Using Mitometer
addpath(strcat(pwd,'/Simulation'))

lengthValues = [1,2,4,8,16];
intValues = [1,2,4,8,9,10];

numTrials = 3;
SNR = zeros(length(lengthValues),length(intValues),numTrials);
meanCompNum = zeros(length(lengthValues),length(intValues),numTrials);
medArea = zeros(length(lengthValues),length(intValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;
    
    noiseLevel = 0;
    for noiseTest = 1:length(intValues)
        lengthNum = 0;
        noiseLevel = noiseLevel+1;
        for lengthTest = 1:length(lengthValues)

            lengthNum = lengthNum+1;
            [mitoIm,mitoSim] = generateIntensityMitochondria(lengthValues(lengthTest),intValues(noiseTest));

            comboMitoIm = zeros(512,512,5);

            tempIm1 = ones(512,512,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
%                 tempIm = tempIm1;
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNR(lengthTest,noiseTest,trial) = signalVal/noiseVal;

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
            medArea(lengthTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNum(lengthTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,lengthTest;meanCompNumTemp,intValues(noiseTest),lengthValues(lengthTest)]
            meanCompNum(:,:,trial)
            trial
        end
    end
end
%%
% Using small medfilt + gaussthresh algorithm

addpath(strcat(pwd,'/Simulation'))

lengthValues = [1,2,4,8,16];
intValues = [1,2,4,8,9,10];

numTrials = 3;
SNRMF3 = zeros(length(lengthValues),length(intValues),numTrials);
meanCompNumMF3 = zeros(length(lengthValues),length(intValues),numTrials);
medAreaMF3 = zeros(length(lengthValues),length(intValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;
    
    noiseLevel = 0;
    for noiseTest = 1:length(intValues)
        lengthNum = 0;
        noiseLevel = noiseLevel+1;
        for lengthTest = 1:length(lengthValues)

            lengthNum = lengthNum+1;
            [mitoIm,mitoSim] = generateIntensityMitochondria(lengthValues(lengthTest),intValues(noiseTest));

            comboMitoIm = zeros(512,512,5);

            tempIm1 = ones(512,512,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
%                 tempIm = tempIm1;
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNRMF3(lengthTest,noiseTest,trial) = signalVal/noiseVal;

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
            medAreaMF3(lengthTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNumMF3(lengthTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,lengthTest;meanCompNumTemp,intValues(noiseTest),lengthValues(lengthTest)]
            meanCompNumMF3(:,:,trial)
            trial
        end
    end
end



% Using big medfilt + gaussthresh algorithm

addpath(strcat(pwd,'/Simulation'))

lengthValues = [1,2,4,8,16];
intValues = [1,2,4,8,9,10];

numTrials = 3;
SNRMF14 = zeros(length(lengthValues),length(intValues),numTrials);
meanCompNumMF14 = zeros(length(lengthValues),length(intValues),numTrials);
medAreaMF14 = zeros(length(lengthValues),length(intValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;
    
    noiseLevel = 0;
    for noiseTest = 1:length(intValues)
        lengthNum = 0;
        noiseLevel = noiseLevel+1;
        for lengthTest = 1:length(lengthValues)

            lengthNum = lengthNum+1;
            [mitoIm,mitoSim] = generateIntensityMitochondria(lengthValues(lengthTest),intValues(noiseTest));

            comboMitoIm = zeros(512,512,5);

            tempIm1 = ones(512,512,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
%                 tempIm = tempIm1;
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNRMF14(lengthTest,noiseTest,trial) = signalVal/noiseVal;

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
            medAreaMF14(lengthTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNumMF14(lengthTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,lengthTest;meanCompNumTemp,intValues(noiseTest),lengthValues(lengthTest)]
            meanCompNumMF14(:,:,trial)
            trial
        end
    end
end


% Using medfilt algorithm + small gauss + thresh algorithm

addpath(strcat(pwd,'/Simulation'))

addpath(strcat(pwd,'/Simulation'))

lengthValues = [1,2,4,8,16];
intValues = [1,2,4,8,9,10];

numTrials = 3;
SNRGF3 = zeros(length(lengthValues),length(intValues),numTrials);
meanCompNumGF3 = zeros(length(lengthValues),length(intValues),numTrials);
medAreaGF3 = zeros(length(lengthValues),length(intValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;
    
    noiseLevel = 0;
    for noiseTest = 1:length(intValues)
        lengthNum = 0;
        noiseLevel = noiseLevel+1;
        for lengthTest = 1:length(lengthValues)

            lengthNum = lengthNum+1;
            [mitoIm,mitoSim] = generateIntensityMitochondria(lengthValues(lengthTest),intValues(noiseTest));

            comboMitoIm = zeros(512,512,5);

            tempIm1 = ones(512,512,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
%                 tempIm = tempIm1;
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNRGF3(lengthTest,noiseTest,trial) = signalVal/noiseVal;

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
            medAreaGF3(lengthTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNumGF3(lengthTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,lengthTest;meanCompNumTemp,intValues(noiseTest),lengthValues(lengthTest)]
            meanCompNumGF3(:,:,trial)
            trial
        end
    end
end

% Using medfilt algorithm + big gauss + thresh algorithm

addpath(strcat(pwd,'/Simulation'))

addpath(strcat(pwd,'/Simulation'))

lengthValues = [1,2,4,8,16];
intValues = [1,2,4,8,9,10];

numTrials = 3;
SNRGF5 = zeros(length(lengthValues),length(intValues),numTrials);
meanCompNumGF5 = zeros(length(lengthValues),length(intValues),numTrials);
medAreaGF5 = zeros(length(lengthValues),length(intValues),numTrials);

for trial = 1:numTrials
    mitoLength = 2;
    vert = 0;
    
    noiseLevel = 0;
    for noiseTest = 1:length(intValues)
        lengthNum = 0;
        noiseLevel = noiseLevel+1;
        for lengthTest = 1:length(lengthValues)

            lengthNum = lengthNum+1;
            [mitoIm,mitoSim] = generateIntensityMitochondria(lengthValues(lengthTest),intValues(noiseTest));

            comboMitoIm = zeros(512,512,5);

            tempIm1 = ones(512,512,'uint8');
            for j = 1:10
                tempIm = imgaussfilt(imnoise(tempIm1,'poisson'),0.3);
%                 tempIm = tempIm1;
                comboMitoImTemp = sum(mitoIm,3);
                comboMitoImTemp = imgaussfilt(uint8(comboMitoImTemp),1);
                comboMitoIm(:,:,j) = comboMitoImTemp+tempIm;
            end

            noiseVal = std(double(tempIm(:)));
            signalVal = mean(comboMitoImTemp(comboMitoImTemp>0));
            SNRtemp = signalVal/noiseVal;
            SNRGF5(lengthTest,noiseTest,trial) = signalVal/noiseVal;

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
            medAreaGF5(lengthTest,noiseTest,trial) = median(areaList);
            meanCompNumTemp = mean(numComps);
            meanCompNumGF5(lengthTest,noiseTest,trial) = median(numComps);

            [medAreaTemp,noiseTest,lengthTest;meanCompNumTemp,intValues(noiseTest),lengthValues(lengthTest)]
            meanCompNumGF5(:,:,trial)
            trial
        end
    end
end