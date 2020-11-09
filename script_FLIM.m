%%
addpath(strcat(pwd,'/FLIM'))

clear
close all

sprintf("Choose your TMRM files in order")
TMRM = Save_R64();

sprintf("Choose your NADH files in the same order")
NADH = Save_R64();
%%
aveNum = 0;

intNADHthresh = 2;
noiseThreshold = 0.25;


TMRMint = zeros(size(TMRM(1).intensity,1),size(TMRM(1).intensity,2),length(TMRM));
NADHint = zeros(size(NADH(1).intensity,1),size(NADH(1).intensity,2),length(NADH));

close all
figure
for i = 1:length(TMRM)
    TMRMint(:,:,i) = uint8(TMRM(i).intensity);
    NADHint(:,:,i) = uint8(NADH(i).intensity);
    
    imagesc(NADHint(:,:,i))
    caxis([0,255])
    pause(0.2)
end
ImAll{1} = TMRMint;
ImNADH{1} = NADHint;

if aveNum
    
    TMRMave = zeros(size(TMRM(1).intensity,1),size(TMRM(1).intensity,2),round(length(TMRM)/aveNum));
    NADHave = zeros(size(NADH(1).intensity,1),size(NADH(1).intensity,2),round(length(NADH)/aveNum));

    % For averaging frames
    for i = 1:size(TMRMint,3)/aveNum
        TMRMave(:,:,i) = mean(TMRMint(:,:,1+aveNum*(i-1):aveNum*i),3);
        NADHave(:,:,i) = mean(NADHint(:,:,1+aveNum*(i-1):aveNum*i),3);
    end
    ImAll{1} = TMRMave;
    ImNADH{1} = NADHave;
    
    close all
    figure
    for i = 2:size(TMRMave,3)
        imagesc(TMRMave(:,:,i))
        caxis([0,255])
        pause(0.05)
    end
end



% close all
% figure
% for i = 2:length(NADH)
% imagesc(ImAll{1}(:,:,i))
% caxis([0,255])
% pause(0.2)
% end

% User Parameters
micronPerPixel = 0.17;
minArea = 0.3;
maxArea = 200;
secondsPerFrame = 5;
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

%
% thr = 10;

%threshold image to remove adjacent mitochondria connections and bg.
ImMask = thresholdImage(ImGaussFiltered,thr);

%apply mask
ImMaskThresholded = areaThreshold(ImMask,minArea/(micronPerPixel^2),maxArea/(micronPerPixel^2));
Im = ImMaskThresholded.*ImOG;

close all
figure
for i = 2:size(Im,3)
imagesc(Im(:,:,i))
caxis([0,255])
pause(0.2)
end

%Preliminary tracking to set weighting
weightsPrelim = ones(1,6)*1/6;
trackPrelim = trackMitochondria(Im,ImOG,weightsPrelim,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds);

%Get the weights
[weights,allPerfectTracks] = getTrackingWeights(trackPrelim);

%Real tracking
[trackAll,mitoCoM,extra] = trackMitochondria(Im,ImOG,weights,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds);
%Get the weights
[weights,track] = getTrackingWeights(trackAll);


% %plot perfect track centroids
% close all
% figure
% imagesc(ImOG(:,:,1))
% colormap gray
% axis image
% hold on
% for i = 1:length(track)
%     if ~isempty(track(i).frame)
%         plotWeightedCentroid(track(i))
% %         plotFissionEvents(allPerfectTracks2,i)
% %         plotFusionEvents(allPerfectTracks2,i)
%     end
% end
% 
% %plot track centroids
% close all
% figure
% imagesc(ImOG(:,:,1))
% colormap gray
% % caxis([-0.05,1]) 
% axis image
% hold on
% for i = 1:length(track)
%     if length(track(i).frame)>0
%         plotWeightedCentroid(track(i))
%         plotFissionEvents(track,i)
% %         plotFusionEvents(track,i)
%     end
% end
% 


%%Mask Flim images

if aveNum
    
    TMRMInttemp = zeros(size(TMRM(1).intensity,1),size(TMRM(1).intensity,2),round(length(TMRM)/aveNum));
    NADHInttemp = zeros(size(NADH(1).intensity,1),size(NADH(1).intensity,2),round(length(NADH)/aveNum));
    TMRMgtemp = zeros(size(TMRM(1).intensity,1),size(TMRM(1).intensity,2),round(length(TMRM)/aveNum));
    NADHgtemp = zeros(size(NADH(1).intensity,1),size(NADH(1).intensity,2),round(length(NADH)/aveNum));
    TMRMstemp = zeros(size(TMRM(1).intensity,1),size(TMRM(1).intensity,2),round(length(TMRM)/aveNum));
    NADHstemp = zeros(size(NADH(1).intensity,1),size(NADH(1).intensity,2),round(length(NADH)/aveNum));

    
    for i = 1:length(TMRM)
        TMRMint(:,:,i) = uint8(TMRM(i).intensity);
        NADHint(:,:,i) = uint8(NADH(i).intensity);
        TMRMgtemp1(:,:,i) = uint8(TMRM(i).g);
        NADHgtemp1(:,:,i) = uint8(NADH(i).g);
        TMRMstemp1(:,:,i) = uint8(TMRM(i).s);
        NADHstemp1(:,:,i) = uint8(NADH(i).s);
        TMRMtptemp1(:,:,i) = uint8(TMRM(i).tauPhase);
    end

    % For averaging frames
    for i = 1:size(TMRMint,3)/aveNum
        TMRMInttemp(:,:,i) = mean(TMRMint(:,:,1+aveNum*(i-1):aveNum*i),3)*aveNum;
        NADHInttemp(:,:,i) = mean(NADHint(:,:,1+aveNum*(i-1):aveNum*i),3)*aveNum;
        TMRMgtemp(:,:,i) = mean(TMRMgtemp1(:,:,1+aveNum*(i-1):aveNum*i),3);
        NADHgtemp(:,:,i) = mean(NADHgtemp1(:,:,1+aveNum*(i-1):aveNum*i),3);
        TMRMstemp(:,:,i) = mean(TMRMstemp1(:,:,1+aveNum*(i-1):aveNum*i),3);
        NADHstemp(:,:,i) = mean(NADHstemp1(:,:,1+aveNum*(i-1):aveNum*i),3);
        TMRMtptemp(:,:,i) = mean(TMRMtptemp1(:,:,1+aveNum*(i-1):aveNum*i),3);
        NADHnew(i).intensity = NADHInttemp(:,:,i);
        TMRMnew(i).intensity = TMRMInttemp(:,:,i);
        NADHnew(i).g = NADHgtemp(:,:,i);
        NADHnew(i).s = NADHstemp(:,:,i);
        TMRMnew(i).g = TMRMgtemp(:,:,i);
        TMRMnew(i).s = TMRMstemp(:,:,i);
        TMRMnew(i).tauPhase = TMRMstemp(:,:,i);
    end
    
    ImAll{1} = TMRMave;
    ImNADH{1} = NADHave;
    NADH = NADHnew;
    TMRM = TMRMnew;
end

[NADHintensity,NADHg,NADHs] = SaveIntGS(NADH,0);
[TMRMintensity,TMRMg,TMRMs] = SaveIntGS(TMRM,0);

TMRMtauP = zeros(size(TMRMintensity));

for i = 1:length(TMRM)
    TMRMtauP(:,:,i) = TMRM(i).tauPhase;
end

GzerosTMRM = TMRMg;
SzerosTMRM = TMRMs;

TMRMg(TMRMg==0) = nan;
TMRMs(TMRMs==0) = nan;

GzerosNADH = NADHg;
SzerosNADH = NADHs;

NADHg(NADHg==0) = nan;
NADHs(NADHs==0) = nan;

% for i = 1:size(NADHg,3)/5
%     NADHintensityAve = mean(NADHintensity(:,:,1+5*(i-1):5*i),3,'omitnan');
%     NADHgAve(:,:,i) = mean(NADHg(:,:,1+5*(i-1):5*i),3,'omitnan');
%     NADHsAve(:,:,i) = mean(NADHs(:,:,1+5*(i-1):5*i),3,'omitnan');
% end
% NADHintensity = NADHintensityAve;
% NADHg = NADHgAve;
% NADHs = NADHsAve;

% NADHs = NADHs*2;
numTracks = length(track);

% for i = 1:size(NADHg,3)
%     NADHg(:,:,i) = medfilt2(NADHg(:,:,i),[5,5]);
%     NADHs(:,:,i) = medfilt2(NADHs(:,:,i),[5,5]);
%     NADHintensity(:,:,i) = medfilt2(NADHintensity(:,:,i),[5,5]);
% end
labelImage = cell(1,numTracks);
labelNADHg = cell(1,numTracks);
labelNADHs = cell(1,numTracks);
labelNADHint = cell(1,numTracks);
labelTMRMg = cell(1,numTracks);
labelTMRMs = cell(1,numTracks);
labelTMRMint = cell(1,numTracks);
labelTMRMtauP = cell(1,numTracks);
fb = cell(1,numTracks);
averageNADHg = zeros(1,numTracks);
averageNADHint = zeros(1,numTracks);
sumNADHint = zeros(1,numTracks);
averageNADHs = zeros(1,numTracks);
stdNADHg = zeros(1,numTracks);
stdNADHs = zeros(1,numTracks);
stdNADHint = zeros(1,numTracks);
averageTMRMg = zeros(1,numTracks);
averageTMRMint = zeros(1,numTracks);
averageTMRMs = zeros(1,numTracks);
averageTMRMtauP = zeros(1,numTracks);
stdTMRMg = zeros(1,numTracks);
stdTMRMs = zeros(1,numTracks);
stdTMRMint = zeros(1,numTracks);
stdTMRMtauP = zeros(1,numTracks);
fFree = zeros(1,numTracks);
fBound = zeros(1,numTracks);
fNoise = zeros(1,numTracks);
fFreeStd = zeros(1,numTracks);
fBoundStd = zeros(1,numTracks);
fNoiseStd = zeros(1,numTracks);

for trackNum = 1:numTracks    
    labelImage{trackNum} = plotLables(track,extra.L,trackNum);
    
    %TMRM
    labelTMRMg{trackNum} = labelImage{trackNum}.*TMRMg;
    labelTMRMs{trackNum} = labelImage{trackNum}.*TMRMs;
    labelTMRMint{trackNum} = labelImage{trackNum}.*TMRMintensity;
    labelTMRMtauP{trackNum} = labelImage{trackNum}.*TMRMtauP;
    averageTMRMtauP(trackNum) = mean(labelTMRMtauP{trackNum}(labelTMRMtauP{trackNum}>0),'all');
    stdTMRMtauP(trackNum) = std(labelTMRMtauP{trackNum}(labelTMRMtauP{trackNum}>0),[],'all');
    averageTMRMint(trackNum) = mean(labelTMRMint{trackNum}(labelTMRMint{trackNum}>0),'all');
    stdTMRMint(trackNum) = std(labelTMRMint{trackNum}(labelTMRMint{trackNum}>0),[],'all');
    averageTMRMg(trackNum) = mean(labelTMRMg{trackNum}(labelTMRMg{trackNum}>0),'all');
    stdTMRMg(trackNum) = std(labelTMRMg{trackNum}(labelTMRMg{trackNum}>0),[],'all');
    averageTMRMs(trackNum) = mean(labelTMRMs{trackNum}(labelTMRMs{trackNum}>0),'all');
    stdTMRMs(trackNum) = std(labelTMRMs{trackNum}(labelTMRMs{trackNum}>0),[],'all');
    
    
    %NADH
    labelNADHg{trackNum} = labelImage{trackNum}.*NADHg;
    labelNADHs{trackNum} = labelImage{trackNum}.*NADHs;
    labelNADHint{trackNum} = labelImage{trackNum}.*NADHintensity;
    averageNADHint(trackNum) = mean(labelNADHint{trackNum}(labelNADHint{trackNum}>0),'all');
    sumNADHint(trackNum) = sum(labelNADHint{trackNum},'all','omitnan');
    stdNADHint(trackNum) = std(labelNADHint{trackNum}(labelNADHint{trackNum}>0),[],'all');
    averageNADHg(trackNum) = mean(labelNADHg{trackNum}(labelNADHg{trackNum}>0),'all');
    stdNADHg(trackNum) = std(labelNADHg{trackNum}(labelNADHg{trackNum}>0),[],'all');
    averageNADHs(trackNum) = mean(labelNADHs{trackNum}(labelNADHs{trackNum}>0),'all');
    stdNADHs(trackNum) = std(labelNADHs{trackNum}(labelNADHs{trackNum}>0),[],'all');
    [frac1,frac2,frac3] = fractionCalc(labelNADHs{trackNum}(labelNADHg{trackNum}>0),labelNADHg{trackNum}(labelNADHg{trackNum}>0),0.4e-09,3.2e-09,80000000);
    frac1 = frac1(abs(frac3)<noiseThreshold);
    frac2 = frac2(abs(frac3)<noiseThreshold);
    frac3 = frac3(abs(frac3)<noiseThreshold);
    fb{trackNum}(:,1) = frac1./(frac1+frac2);
    fb{trackNum}(:,2) = frac2./(frac1+frac2);
    fb{trackNum}(:,3) = frac3;
    fFree(trackNum) = nanmean(frac1./(frac1+frac2));
    fFreeStd(trackNum) = nanstd(frac1./(frac1+frac2));
    fBound(trackNum) = nanmean(frac2./(frac1+frac2));
    fBoundStd(trackNum) = nanstd(frac2./(frac1+frac2));
    fNoise(trackNum) = nanmean(frac3);
    fNoiseStd(trackNum) = nanstd(frac3);

end


%cytoplasm,mito

allg = NADHg(NADHintensity>intNADHthresh);
alls = NADHs(NADHintensity>intNADHthresh);

[allfrac1,allfrac2,allfrac3] = fractionCalc(alls,allg,0.4e-09,3.2e-09,80000000);

allfFree = allfrac1./(allfrac1+allfrac2);
allfFree(abs(allfrac3)>noiseThreshold) = nan;
allfFree = allfFree(~isnan(allfFree));
allfBound = allfrac2./(allfrac1+allfrac2);
allfBound(abs(allfrac3)>noiseThreshold) = nan;
allfBound = allfBound(~isnan(allfBound));

allfFreemean = nanmean(allfFree);
allfBoundmean = nanmean(allfBound);

allfFreestd = nanstd(allfFree);
allfBoundstd = nanstd(allfBound);

cytoNADHg = ~ImMask.*NADHg;
mitoNADHg = ImMaskThresholded.*NADHg;
cytoNADHs = ~ImMask.*NADHs;
mitoNADHs = ImMaskThresholded.*NADHs;

cytoNADHinttemp = ~ImMask.*NADHintensity;
cytoNADHint = cytoNADHinttemp;
cytoNADHint(cytoNADHinttemp<intNADHthresh) = nan;
cytoNADHg(cytoNADHinttemp<intNADHthresh) = nan;
cytoNADHs(cytoNADHinttemp<intNADHthresh) = nan;

mitoNADHinttemp = ImMaskThresholded.*NADHintensity;
mitoNADHint = mitoNADHinttemp;
mitoNADHint(mitoNADHinttemp<intNADHthresh) = nan;
mitoNADHg(mitoNADHinttemp<intNADHthresh) = nan;
mitoNADHs(mitoNADHinttemp<intNADHthresh) = nan;

averageCytoNADHint = mean(cytoNADHint(cytoNADHint>0),'all');
averageCytoNADHg = mean(cytoNADHg(cytoNADHg>0),'all','omitnan');
averageCytoNADHs = mean(cytoNADHs(cytoNADHs>0),'all','omitnan');

averageMitoNADHint = mean(mitoNADHint(mitoNADHint>0),'all');
averageMitoNADHg = mean(mitoNADHg(mitoNADHg>0),'all','omitnan');
averageMitoNADHs = mean(mitoNADHs(mitoNADHs>0),'all','omitnan');

[cytofrac1,cytofrac2,cytofrac3] = fractionCalc(cytoNADHs(~isnan(cytoNADHs)),cytoNADHg(~isnan(cytoNADHg)),0.4e-09,3.2e-09,80000000);
[mitofrac1,mitofrac2,mitofrac3] = fractionCalc(mitoNADHs(~isnan(mitoNADHs)),mitoNADHg(~isnan(mitoNADHg)),0.4e-09,3.2e-09,80000000);

cytofFree = cytofrac1./(cytofrac1+cytofrac2);
cytofFree(abs(cytofrac3)>noiseThreshold) = nan;
cytofFree = cytofFree(~isnan(cytofFree));
cytofBound = cytofrac2./(cytofrac1+cytofrac2);
cytofBound(abs(cytofrac3)>noiseThreshold) = nan;
cytofBound = cytofBound(~isnan(cytofBound));

mitofFree = mitofrac1./(mitofrac1+mitofrac2);
mitofFree(abs(mitofrac3)>noiseThreshold) = nan;
mitofFree = mitofFree(~isnan(mitofFree));
mitofBound = mitofrac2./(mitofrac1+mitofrac2);
mitofBound(abs(mitofrac3)>noiseThreshold) = nan;
mitofBound = mitofBound(~isnan(mitofBound));

cytofFreemean = nanmean(cytofFree);
mitofFreemean = nanmean(mitofFree);
cytofBoundmean = nanmean(cytofBound);
mitofBoundmean = nanmean(mitofBound);

cytofFreestd = nanstd(cytofFree);
mitofFreestd = nanstd(mitofFree);
cytofBoundstd = nanstd(cytofBound);
mitofBoundstd = nanstd(mitofBound);

cytofFreenum = numel(cytofFree(~isnan(cytofFree)));
mitofFreenum = numel(mitofFree(~isnan(mitofFree)));
cytofBoundnum = numel(cytofBound(~isnan(cytofBound)));
mitofBoundnum = numel(mitofBound(~isnan(mitofBound)));

cytofNoisestd = nanstd(cytofrac3);
mitofNoisestd = nanstd(mitofrac3);
%cytomitoplasm done

GScoordsNADH = [averageNADHg',averageNADHs'];
GSstdcoordsNADH = [stdNADHg',stdNADHs'];

GScoordsTMRM = [averageTMRMg',averageTMRMs'];
GSstdcoordsTMRM = [stdTMRMg',stdTMRMs'];

script_FLIManalysis

%% Mean Intensity Cutoff NADH
close all

IntCutoff = -1;

% trackToView = (1:length(track));
trackToView = 17;
numTracks = length(trackToView);
meanInt = zeros(1,numTracks);

% for trackNum = 1:numTracks
%     meanInt(trackNum) = mean(track(trackNum).MeanIntensity);
% end
for trackNum = ThrTracks
    meanInt(trackNum) = mean(track(trackNum).MeanIntensity);
end

% GScutoff = [GScoords(meanInt>IntCutoff,1),GScoords(meanInt>IntCutoff,2)];
% GSstdcutoff = [GSstdcoords(meanInt>IntCutoff,1),GSstdcoords(meanInt>IntCutoff,2)];
GScutoff = [GScoordsNADH(trackToView,1),GScoordsNADH(trackToView,2)];
GSstdcutoff = [GSstdcoordsNADH(trackToView,1),GSstdcoordsNADH(trackToView,2)];


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

medfiltMeanG = medfilt2(mean(GzerosNADH,3),[5,5]);
medfiltMeanS = medfilt2(mean(SzerosNADH,3),[5,5]);
meanInt = mean(NADHintensity,3);
intThresh = 0;
threshG = medfiltMeanG(meanInt>=intThresh);
threshS = medfiltMeanS(meanInt>=intThresh);

histogram2(threshG,threshS,[round(numel(threshG)/50),round(numel(threshS)/50)],'FaceColor','flat')
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


%% Mean Intensity Cutoff TMRM
close all

IntCutoff = -1;
% trackToView = (1:length(track));
trackToView = ThrTracks;
numTracks = length(trackToView);
meanInt = zeros(1,numTracks);

% for trackNum = 1:numTracks
%     meanInt(trackNum) = mean(track(trackNum).MeanIntensity);
% end
for trackNum = ThrTracks
    meanInt(trackNum) = mean(track(trackNum).MeanIntensity);
end


% GScutoff = [GScoords(meanInt>IntCutoff,1),GScoords(meanInt>IntCutoff,2)];
% GSstdcutoff = [GSstdcoords(meanInt>IntCutoff,1),GSstdcoords(meanInt>IntCutoff,2)];
GScutoff = [GScoordsTMRM(trackToView,1),GScoordsTMRM(trackToView,2)];
GSstdcutoff = [GSstdcoordsTMRM(trackToView,1),GSstdcoordsTMRM(trackToView,2)];


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
medfiltMeanG = medfilt2(mean(GzerosTMRM,3),[5,5]);
medfiltMeanS = medfilt2(mean(SzerosTMRM,3),[5,5]);
meanInt = mean(TMRMintensity,3);
intThresh = 10;
threshG = medfiltMeanG(meanInt>=intThresh);
threshS = medfiltMeanS(meanInt>=intThresh);

histogram2(threshG,threshS,[round(numel(threshG)/50),round(numel(threshS)/50)],'FaceColor','flat')
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

% trackToView = 1:length(ThrTracks);
addedAll = zeros(size(TMRM(1).intensity(:,:,1),1),size(TMRM(1).intensity(:,:,1),2));

trackToView = ThrTracks;
trackToView = 17;

for trackNum = trackToView
    labelImage = plotLables(track,extra.L,trackNum);
    for frameNum = 1:size(labelImage,3)
        currIm = labelImage(:,:,frameNum);
        addedAll(currIm>0) = frameNum;
    end
end

imagesc(addedAll)
axis image off
colormap jet


%%
%plot track centroids
close all
figure
imagesc(ImOG(:,:,1))
colormap gray
% caxis([-0.05,1]) 
axis image
hold on
for i = 1:length(track)
    if length(track(i).frame)>3
        plotWeightedCentroid(track(i))
        plotFissionEvents(track,i)
        plotFusionEvents(track,i)
    end
end
