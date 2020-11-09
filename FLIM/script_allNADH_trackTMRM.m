
addpath(strcat(pwd,'/FLIM'))

clear
close all

sprintf("Choose your files in order")
[Im,Names,Path] = Open_R64();
%%
clear TMRM
clear NADH

NADHind = 1;
TMRMind = 1;

NADH = struct();
TMRM = struct();

omega = 80000000;

for imageNum = 1:size(Im,2) 
    
    %cut off the edges because they are weird.
    im(imageNum).intensity = Im{1,imageNum}(2:end-2,2:end-2,1);
    if isempty(min(im(imageNum).intensity(im(imageNum).intensity>0)))
        continue
    end
    im(imageNum).intensity = im(imageNum).intensity/min(im(imageNum).intensity(im(imageNum).intensity>0));
    
    
    im(imageNum).phase = Im{1,imageNum}(2:end-2,2:end-2,2);
    im(imageNum).mod = Im{1,imageNum}(2:end-2,2:end-2,3);
    im(imageNum).phase2 = Im{1,imageNum}(2:end-2,2:end-2,2);
    im(imageNum).mod2 = Im{1,imageNum}(2:end-2,2:end-2,3);
    
    im(imageNum).tauPhase = tand(Im{1,imageNum}(2:end-2,2:end-2,2))/(omega);
    im(imageNum).tauMod = real(sqrt(((1./(Im{1,imageNum}(2:end-2,2:end-2,3))).^2)-1)./(omega));
    
    im(imageNum).g = Im{1,imageNum}(2:end-2,2:end-2,3).*cosd(Im{1,imageNum}(2:end-2,2:end-2,2));
    im(imageNum).s = Im{1,imageNum}(2:end-2,2:end-2,3).*sind(Im{1,imageNum}(2:end-2,2:end-2,2));
    
    im(imageNum).g2 = Im{1,imageNum}(2:end-2,2:end-2,5).*cosd(Im{1,imageNum}(2:end-2,2:end-2,4));
    im(imageNum).s2 = Im{1,imageNum}(2:end-2,2:end-2,5).*sind(Im{1,imageNum}(2:end-2,2:end-2,4));
    
    im(imageNum).name = Names{imageNum}; 
    im(imageNum).path = Path;
    
    if strcmp(Names{imageNum}(end-10),'1')
        if NADHind == 1
            NADH = im(imageNum);
            NADHind = NADHind+1;
        else
            NADH(NADHind) = im(imageNum);
            NADHind = NADHind+1;
        end
    elseif strcmp(Names{imageNum}(end-10),'2')
        if TMRMind == 1
            TMRM = im(imageNum);
            TMRMind = TMRMind+1;
        else
            TMRM(TMRMind) = im(imageNum);
            TMRMind = TMRMind+1;
        end
    end
    
end
%%
NADH = Save_R64();
sprintf("Choose your TMRM files in order")
TMRM = Save_R64();

%%
close all
figure
for i = 1:length(TMRM)
    imagesc(TMRM(i).intensity)
    axis image
%     caxis([0,40])
    pause(0.05)
end

%%
TMRMint = zeros(size(TMRM(1).intensity,1),size(TMRM(1).intensity,2),length(TMRM));
NADHint = NADH.intensity;

for i = 1:length(TMRM)
    TMRMint(:,:,i) = uint8(TMRM(i).intensity);
end

% User Parameters
micronPerPixel = 0.18;
minArea = 0.3;
maxArea = 200;
secondsPerFrame = 2.52;
distThreshMicrons = 1;
frameThreshSeconds = 15;

% User Image
ImOG = TMRMint;

%remove diffuse bg and normalize image intensity
[ImBgRemoved,ImMinMed,ImMedFiltered] = diffuseBgRemove(ImOG,(minArea)/micronPerPixel,(1)/micronPerPixel);

%find optimal gaussian filter std Dev and hard threshold value
[sig,thr,costs] = optimizeSigmaThresh(ImBgRemoved);

%gaussian filter for low-pass filter and smooth discontinuities.
ImGaussFiltered = gaussFilter(ImBgRemoved,sig);
%%
thr = 1;
%%
%threshold image to remove adjacent mitochondria connections and bg.
ImMask = thresholdImage(ImGaussFiltered,thr);

%apply mask
ImMaskThresholded = areaThreshold(ImMask,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));
Im = ImMaskThresholded.*ImOG;

close all
figure
for i = 1:size(TMRMint,3)
    imagesc(Im(:,:,i))
    axis image
%     caxis([0,40])
    pause(0.02)
end
%%
%Preliminary tracking to set weighting
weightsPrelim = ones(1,6)*1/6;
trackPrelim = trackMitochondria(Im,ImOG,weightsPrelim,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds,1,1);

%Get the weights
[weights,allPerfectTracks] = getTrackingWeights(trackPrelim);

%Real tracking
[trackAll,mitoCoM,extra] = trackMitochondria(Im,ImOG,weights,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds,1,1);
%Get the weights
[weights,track] = getTrackingWeights(trackAll);

%%
%plot perfect track centroids
close all
figure
imagesc(ImOG(:,:,1))
colormap gray
axis image
hold on
for i = 1:length(track)
    if ~isempty(track(i).frame)
        plotWeightedCentroid(track(i))
%         plotFissionEvents(track,i)
%         plotFusionEvents(track,i)
    end
end

%%
%plot track centroids
close all
figure
imagesc(ImOG(:,:,1))
colormap gray
% caxis([-0.05,1]) 
axis image
hold on
for i = 1:length(trackAll)
    if length(trackAll(i).frame)>3
        plotWeightedCentroid(trackAll(i))
        plotFissionEvents(trackAll,i)
        plotFusionEvents(trackAll,i)
    end
end

%% NADH phasor for individual mitochondria
numTracks = length(track);
sumGWeightedNormalized = zeros(1,numTracks);
sumSWeightedNormalized = zeros(1,numTracks);
averageGAll = cell(1,numTracks);
averageSAll = cell(1,numTracks);
averageG = zeros(1,numTracks);
averageS = zeros(1,numTracks);
stdG = zeros(1,numTracks);
stdS = zeros(1,numTracks);

for trackNum = 1:numTracks
    totalPixelReps = zeros(size(NADH.g,1)*size(NADH.g,2),1);
    for indexNum = 1:length(track(trackNum).Area)
        if ~iscell(track(trackNum).PixelIdxList)
            track(trackNum).PixelIdxList = {track(trackNum).PixelIdxList};
        end
        totalPixelReps(track(trackNum).PixelIdxList{indexNum}) = totalPixelReps(track(trackNum).PixelIdxList{indexNum})+1;
        averageGAll{trackNum}(indexNum) = nanmean(NADH.g(track(trackNum).PixelIdxList{indexNum}));
        averageSAll{trackNum}(indexNum) = nanmean(NADH.s(track(trackNum).PixelIdxList{indexNum}));
    end
    gWeighted = NADH.g(:).*totalPixelReps;
    sWeighted = NADH.s(:).*totalPixelReps;
    sumGWeightedNormalized(trackNum) = sum(gWeighted)/sum(totalPixelReps);
    sumSWeightedNormalized(trackNum) = sum(sWeighted)/sum(totalPixelReps);
    
    averageG(trackNum) = nanmean(averageGAll{trackNum});
    stdG(trackNum) = nanstd(averageGAll{trackNum});
    averageS(trackNum) = nanmean(averageSAll{trackNum});
    stdS(trackNum) = nanstd(averageSAll{trackNum});

end
close all
figure

%Plot g and s
%plot g and s of tracks above mean intensity cutoff
eb(1) = errorbar(averageG,averageS,stdG, 'horizontal', 'LineStyle', 'none');
hold on
eb(2) = errorbar(averageG,averageS,stdS, 'vertical', 'LineStyle', 'none');
scatter(averageG,averageS,10,'r','filled','MarkerEdgeColor','k','LineWidth',1)
set(eb, 'color', 'k', 'LineWidth', 0.5,'CapSize',1)
colormap jet
grid on
view(2)
hold on
% scatter(sumGWeightedNormalized,sumSWeightedNormalized)

%Draw semi-circle
th = linspace(0, pi, 100);
R = 0.5;
x = R*cos(th) + 0.5;
y = R*sin(th);
plot(x,y,'k','LineWidth',2); 
axis equal;
xlim([0,1])
ylim([0,1])


%% TMRM phasor for individual mitochondria pt1

sprintf("Choose your TMRM averaged File")
TMRMAll = Save_R64();
%% pt2
numTracks = length(track);
sumGWeightedNormalizedTMRM = zeros(1,numTracks);
sumSWeightedNormalizedTMRM = zeros(1,numTracks);
averageGAllTMRM = cell(1,numTracks);
averageSAllTMRM = cell(1,numTracks);
averageTpAllTMRM = cell(1,numTracks);
averageGTMRM = zeros(1,numTracks);
averageSTMRM = zeros(1,numTracks);
averageTpTMRM = zeros(1,numTracks);
stdTpTMRM = zeros(1,numTracks);
stdGTMRM = zeros(1,numTracks);
stdSTMRM = zeros(1,numTracks);

for trackNum = 1:numTracks
    totalPixelReps = zeros(size(TMRMAll.g,1)*size(TMRMAll.g,2),1);
    for indexNum = 1:length(track(trackNum).Area)
        totalPixelReps(track(trackNum).PixelIdxList{indexNum}) = totalPixelReps(track(trackNum).PixelIdxList{indexNum})+1;
        averageGAllTMRM{trackNum}(indexNum) = nanmean(TMRMAll.g(track(trackNum).PixelIdxList{indexNum}));
        averageSAllTMRM{trackNum}(indexNum) = nanmean(TMRMAll.s(track(trackNum).PixelIdxList{indexNum}));
        averageTpAllTMRM{trackNum}(indexNum) = nanmean(TMRMAll.tauPhase(track(trackNum).PixelIdxList{indexNum}));
    end
    gWeighted = TMRMAll.g(:).*totalPixelReps;
    sWeighted = TMRMAll.s(:).*totalPixelReps;
    sumGWeightedNormalizedTMRM(trackNum) = sum(gWeighted)/sum(totalPixelReps);
    sumSWeightedNormalizedTMRM(trackNum) = sum(sWeighted)/sum(totalPixelReps);
    
    averageTpTMRM(trackNum) = nanmean(averageTpAllTMRM{trackNum});
    stdTpTMRM(trackNum) = nanstd(averageTpAllTMRM{trackNum});
    averageGTMRM(trackNum) = nanmean(averageGAllTMRM{trackNum});
    stdGTMRM(trackNum) = nanstd(averageGAllTMRM{trackNum});
    averageSTMRM(trackNum) = nanmean(averageSAllTMRM{trackNum});
    stdSTMRM(trackNum) = nanstd(averageSAllTMRM{trackNum});
end

close all
figure

%Plot g and s
%plot g and s of tracks above mean intensity cutoff
eb(1) = errorbar(averageGTMRM,averageSTMRM,stdGTMRM, 'horizontal', 'LineStyle', 'none');
hold on
eb(2) = errorbar(averageGTMRM,averageSTMRM,stdSTMRM, 'vertical', 'LineStyle', 'none');
scatter(averageGTMRM,averageSTMRM,10,'r','filled','MarkerEdgeColor','k','LineWidth',1)
set(eb, 'color', 'k', 'LineWidth', 0.5,'CapSize',1)
colormap jet
grid on
view(2)
hold on

%Draw semi-circle
th = linspace(0, pi, 100);
R = 0.5;
x = R*cos(th) + 0.5;
y = R*sin(th);
plot(x,y,'k','LineWidth',2); 
axis equal;
xlim([0,1])
ylim([0,1])

%% fraction calc

noiseThreshold = 0.5;

[allfrac1,allfrac2,allfrac3] = fractionCalc(averageS,averageG,0.4e-09,3.2e-09,80000000);
% [allfrac1std1,allfrac2std1,allfrac3std1] = fractionCalc(averageS-stdS,averageG+stdG,0.4e-09,3.2e-09,80000000);
% [allfrac1std2,allfrac2std2,allfrac3std2] = fractionCalc(averageS+stdS,averageG-stdG,0.4e-09,3.2e-09,80000000);

allfFree = allfrac1./(allfrac1+allfrac2);
allfFree(abs(allfrac3)>noiseThreshold) = nan;
allfBound = allfrac2./(allfrac1+allfrac2);
allfBound(abs(allfrac3)>noiseThreshold) = nan;

%% Color based on fraction bound

numFrames = length(extra.L);
close all
figure
fractionImAll = zeros(253,253,numFrames);
for frameNum = 1:numFrames
    newImLabel = double(extra.L{frameNum});
    fractionIm = zeros(size(newImLabel,1),size(newImLabel,2));

    for j = 1:length(track)
        if ~sum(track(j).frame == frameNum) || isnan(allfBound(j))
            continue
        end
        frameIdx = find(track(j).frame==frameNum);
        labelIdx = find(newImLabel==track(j).label(frameIdx));
        fractionIm(labelIdx) = allfBound(j);
    end

    imbg = (ImOG(:,:,frameNum)/max(ImOG(:)));
    imbg(fractionIm ~= 0) = 0;

    imagesc(fractionIm)% + imbg)
    fractionImAll(:,:,frameNum) = fractionIm;
    axis image off
    colormap jet

    caxis([0,1])
    pause(0.02)
end

%% Color based on TMRM tauPhase

numFrames = length(extra.L);
close all
figure
tauPhaseImAll = zeros(253,253,numFrames);
for frameNum = 1:numFrames
    newImLabel = double(extra.L{frameNum});
    tauPhaseIm = zeros(size(newImLabel,1),size(newImLabel,2));

    for j = 1:length(track)
        if ~sum(track(j).frame == frameNum) || isnan(allfBound(j))
            continue
        end
        frameIdx = find(track(j).frame==frameNum);
        labelIdx = find(newImLabel==track(j).label(frameIdx));
        tauPhaseIm(labelIdx) = averageTpTMRM(j);
    end

    imbg = (ImOG(:,:,frameNum)/max(ImOG(:)))*max(averageTpTMRM(:));
    imbg(tauPhaseIm ~= 0) = 0;

%     
%     imagesc(tauPhaseIm)
    tauPhaseTemp = tauPhaseIm;
    tauPhaseTemp(tauPhaseTemp>2.5E-9) = 2.5E-9;
    tauPhaseTemp(tauPhaseTemp<1.7E-9) = 1.7E-9;
    tauPhaseIm = tauPhaseTemp - min(tauPhaseTemp(:));
    tauPhaseIm = tauPhaseIm ./ max(tauPhaseIm(:));
    tauPhaseImAll(:,:,frameNum) = tauPhaseIm;
        
    imagesc(tauPhaseImAll(:,:,frameNum))

    axis image off
    colormap jet

%     caxis([1.7E-9,2.5E-9])
    pause(0.02)
end


%% Get dynamics

[dyn,CoMdyn] = getDynamics(0.18,2.52,track,mitoCoM);
for i = 1:length(track)
    allArea(i) = mean(track(i).Area);
    allMajorAxisLength(i) = mean(track(i).MajorAxisLength);
    allEccentricity(i) = mean(track(i).Eccentricity);
    allPerimeter(i) = mean(track(i).Perimeter);
    allSpeed(i) = nanmean(dyn.speed{i});
    sumDistance(i) = sum(dyn.distance{i});
    finalDisplacement(i) = dyn.displacement{i}(end);
    dRatio(i) = finalDisplacement(i)/sumDistance(i);
    allTMRMint(i) = mean(track(i).MeanIntensity);
end
%%
[averageG;stdG;averageS;stdS;averageGTMRM;stdGTMRM;averageSTMRM;stdSTMRM;allfBound;allArea;allMajorAxisLength;allEccentricity;allPerimeter;allSpeed;sumDistance;finalDisplacement;dRatio;allTMRMint]