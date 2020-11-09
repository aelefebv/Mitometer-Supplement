clear
close all
clc

[fileName,path] = uigetfile('*.tif','Select one or more thresholded files','MultiSelect','on');

if ~iscell(fileName)
    fileName = {fileName};
end

numIms = length(fileName);

% [ImAll,~,fileName,path] = openTifs();
% %%

% User Parameters
micronPerPixel = 0.0879;
minArea = 0.3;
maxArea = 200;
secondsPerFrame = 1;
distThreshMicrons = 1;

for imNum = 1:numIms
    clear holdIms
    [ImAll,~,~,~] = openTifs(1,fileName(imNum),path);
      
    strcat("Image number ",num2str(imNum)," of ",num2str(numIms))
    
    % User Image
    ImOG = ImAll{1}(:,:,1:120);
    holdIms(imNum).ImOG = ImOG;
    holdIms(imNum).fileName = fileName{imNum};
    holdIms(imNum).path = path;
    
    %remove diffuse bg and normalize image intensity
    [ImBgRemoved,ImMinMed,ImMedFiltered] = diffuseBgRemove(ImOG,(minArea)/micronPerPixel,(1)/micronPerPixel);

    %find optimal gaussian filter std Dev and hard threshold value
    [sig,thr,costs] = optimizeSigmaThresh(ImBgRemoved,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));
    holdIms(imNum).sig = sig;
    holdIms(imNum).thr = thr;

    %gaussian filter for low-pass filter and smooth discontinuities.
    ImGaussFiltered = gaussFilter(ImBgRemoved,sig);

    %threshold image to remove adjacent mitochondria connections and bg.
    ImMask = thresholdImage(ImGaussFiltered,thr);

    %apply mask
    ImMaskThresholded = areaThreshold(ImMask,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));
    Im = ImMaskThresholded.*ImOG;
    holdIms(imNum).Im = Im;
    
%     Im = tempIm(:,:,1:4);
%     ImOG = tempImOG(:,:,1:4);
    
%     outputFileName = strcat(path,"/segmented_",fileName{imNum});
%     for iii=1:size(Im,3)
%         for iv = 1:size(Im,4)
%             imwrite(Im(:,:,iii,iv), outputFileName, 'WriteMode', 'append',  'Compression','none');
%         end
%     end

    searchTimeVec = 1:5:120;
    for searchTime = 1:length(searchTimeVec)
        frameThreshSeconds = searchTimeVec(searchTime);

        %Preliminary tracking to set weighting
        weightsPrelim = ones(1,6)*1/6;
        trackPrelim = trackMitochondria(Im,ImOG,weightsPrelim,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds,0,1);
        %Get the weights
        [weights,allPerfectTracks] = getTrackingWeights(trackPrelim);

        %Real tracking
        [track,mitoCoM,extra] = trackMitochondria(Im,ImOG,weights,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds,0,1);
        %Get the weights
        [weights2,allPerfectTracks2] = getTrackingWeights(track);

        aveMito = zeros(1,length(extra.mito));
        for i = 1:length(extra.mito)
            aveMito(i) = length(extra.mito{i});
        end

        aveAll = mean(aveMito);

        num = 1;
        big = 0;
        small = 0;
        found = 0;
        test = trackLengthThreshold(track,num);
        while found == 0
            if length(test)>aveAll && small == 0
                big = 1;
                small = 0;
                num = num + 1;
                test = trackLengthThreshold(track,num);
            elseif length(test)<aveAll && big == 0
                small = 1;
                big = 0;
                num = num-1;
                test = trackLengthThreshold(track,num);
            else
                found =1;
            end
        end
        optimalThresh = num;

        dynamics = getDynamics(micronPerPixel,secondsPerFrame,track);

%         holdIms(imNum).optimalThresh = optimalThresh;
%         holdIms(imNum).weights = weights;
%         holdIms(imNum).allPerfectTracks = allPerfectTracks;
%         holdIms(imNum).trackPrelim = trackPrelim;
%         holdIms(imNum).track = track;
%         holdIms(imNum).dynamics = dynamics;
%         holdIms(imNum).mitoCoM = mitoCoM;
%         holdIms(imNum).extra = extra;
%         holdIms(imNum).weights2 = weights2;
%         holdIms(imNum).allPerfectTracks2 = allPerfectTracks2;

%         saveTemp.sig = holdIms(imNum).sig;
%         saveTemp.thr = holdIms(imNum).thr;
%         saveTemp.weights = holdIms(imNum).weights;
%         saveTemp.optimalThresh = holdIms(imNum).optimalThresh;
%         saveTemp.mitoCoM = holdIms(imNum).mitoCoM;
%         saveTemp.track = holdIms(imNum).track;
        
        searchTimeSave(imNum).searchTimeVal(searchTime).track = track;
        imNum
        searchTime
    end
    
    
% %     saveTemp = holdIms(imNum);  
%     saveLoc = strcat(path,'workspace_',fileName{imNum},'.mat');
%     if ~iscell(saveLoc)
%         saveLoc = {saveLoc};
%     end
%     save(saveLoc{1},'saveTemp');
    
end
