%%
addpath(strcat(pwd,'/FLIM'))

clear ImAll
TMRM = Save_R64();
for i = 1:length(TMRM)
TMRMint(:,:,i) = TMRM(i).intensity;
end
ImAll{1} = im2uint8(TMRMint);

%%
close all
figure
for i = 2:length(TMRM)
imagesc(ImAll{1}(:,:,i))
caxis([0,255])
pause(0.02)
end

%%
% User Parameters
micronPerPixel = 0.18;
minArea = 0.3;
maxArea = 1e4;
secondsPerFrame = 3.85;
distThreshMicrons = 1;
frameThreshSeconds = 15;

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
Im = ImMaskThresholded.*ImOG;
%%
%Preliminary tracking to set weighting
weightsPrelim = ones(1,6)*1/6;
trackPrelim = trackMitochondria(Im,ImOG,weightsPrelim,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds);

%Get the weights
[weights,allPerfectTracks] = getTrackingWeights(trackPrelim);

%Real tracking
[track,mitoCoM,extra] = trackMitochondria(Im,ImOG,weights,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds);
%Get the weights
[weights,allPerfectTracks2] = getTrackingWeights(trackPrelim);

%%
%plot perfect track centroids
close all
figure
imagesc(255-ImOG(:,:,1))
colormap gray
axis image
hold on
for i = 1:length(allPerfectTracks2)
    if length(allPerfectTracks2(i).frame)>20
        plotWeightedCentroid(allPerfectTracks2(i))
%         plotFissionEvents(allPerfectTracks2,i)
%         plotFusionEvents(allPerfectTracks2,i)
    end
end
%%
%plot track centroids
close all
figure
imagesc(255-ImOG(:,:,1))
colormap gray
% caxis([-0.05,1]) 
axis image
hold on
for i = 1:length(track)
    if length(track(i).frame)>10
        plotWeightedCentroid(track(i))
        plotFissionEvents(track,i)
        plotFusionEvents(track,i)
    end
end

%%
ImOG = openTifs();
%%
%plot perfect track centroids
close all
figure
imagesc(255-ImOG{1}(:,:,1))
colormap gray
% caxis([0,255])
axis image
hold on
numTrackAll = 0;
[~,prekept] = getTrackingWeights(saveTemp.track);
kept = trackLengthThreshold(saveTemp.track,saveTemp.optimalThresh);

for i = 1:length(kept)
    if length(kept(i).frame)>0
        numTrackAll = numTrackAll+1;
        plotWeightedCentroid(kept(i))
%         i
%         pause(0.2)
        plotFissionEvents(kept,i)
        plotFusionEvents(kept,i)
    end
end
numTrackAll

%% Mask Flim images

NADH = Save_R64();
tempTrack = track;
track = allPerfectTracks2;
[NADHintensity,NADHg,NADHs] = SaveIntGS(NADH,0);
numTracks = length(track);

labelImage = cell(1,numTracks);
labelG = cell(1,numTracks);
labelS = cell(1,numTracks);
averageG = zeros(1,numTracks);
averageS = zeros(1,numTracks);
stdG = zeros(1,numTracks);
stdS = zeros(1,numTracks);

for trackNum = 1:numTracks
    labelImage{trackNum} = plotLables(track,extra.L,trackNum);
    labelG{trackNum} = labelImage{trackNum}.*NADHg;
    labelS{trackNum} = labelImage{trackNum}.*NADHs;
    averageG(trackNum) = mean(labelG{trackNum}(labelG{trackNum}>0 && labelG{trackNum}<1),'all');
    stdG(trackNum) = std(labelG{trackNum}(labelG{trackNum}>0 && labelG{trackNum}<1),[],'all');
    averageS(trackNum) = mean(labelS{trackNum}(labelS{trackNum}>0 && labelS{trackNum}<0.5),'all');
    stdS(trackNum) = std(labelS{trackNum}(labelS{trackNum}>0 && labelS{trackNum}<0.5),[],'all');
end

GScoords = [averageG',averageS'];
GSstdcoords = [stdG',stdS'];

%% Mean Intensity Cutoff
close all

IntCutoff = 1;

meanInt = zeros(1,numTracks);

for trackNum = 1:numTracks
    meanInt(trackNum) = mean(track(trackNum).MeanIntensity);
end

GScutoff = [GScoords(meanInt>IntCutoff,1),GScoords(meanInt>IntCutoff,2)];
GSstdcutoff = [GSstdcoords(meanInt>IntCutoff,1),GSstdcoords(meanInt>IntCutoff,2)];

%plot g and s of tracks above mean intensity cutoff
eb(1) = errorbar(GScutoff(:,1),GScutoff(:,2),GSstdcutoff(:,1), 'horizontal', 'LineStyle', 'none');
hold on
eb(2) = errorbar(GScutoff(:,1),GScutoff(:,2),GSstdcutoff(:,2), 'vertical', 'LineStyle', 'none');
scatter(GScutoff(:,1),GScutoff(:,2),10,'r','filled','MarkerEdgeColor','k','LineWidth',1)

set(eb, 'color', 'k', 'LineWidth', 0.5,'CapSize',1)
colormap jet
grid on


%Draw semi-circle
th = linspace(0, pi, 100);
R = 0.5;
x = R*cos(th) + 0.5;
y = R*sin(th);
plot(x,y,'k','LineWidth',2); 
axis equal;
xlim([0,1])
ylim([0,1])

%%Plot all G and S original
figure

%Plot G and S
histogram2(NADHg(NADHintensity>IntCutoff),NADHs(NADHintensity>IntCutoff),[600,600],'FaceColor','flat')
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
hold on
xlim([0,1])
ylim([0,1])

%% View all tracks
close all
figure

% IntCutoff = 40;

% trackToView = find(meanInt>IntCutoff); %above a certain intensity

addedAll = zeros(size(labelImage,1),size(labelImage,2));

for trackNum = trackToView
    labelImage = plotLables(track,extra.L,trackNum);
    for frameNum = 1:size(labelImage,3)
        currIm = labelImage(:,:,frameNum);
        addedAll(currIm>0) = frameNum;
    end
end

imagesc(addedAll)
colormap jet