close all
clear

addpath(strcat(pwd,'/Simulation'))

densityTrials = [20,40,60,80];
numDensity = length(densityTrials);

comboStats = cell(1,numDensity);
fusionStatsCell = cell(1,numDensity);

for densityNum = 1:numDensity
    numFrames = 20;
    numMitoOG = densityTrials(densityNum);
    frameSize = 512;
    micronPerPixel = 0.138;
    secondsPerFrame = 5;
    lengthMultiplier = 1;
    intensityMultiplier = 1;
    speedMultiplier = 1;
    percentNewMito = 0;
    lostChance = 0.05;
    fusionChance = 1;
    fissionChance = 0;
    noiseMultiplier = 1;

    allMito = cell(1,numFrames);
    mitoIm = cell(1,numFrames);
    newMitoIm = cell(1,numFrames);
    mitoSim = cell(1,numFrames);

    close all
    [mitoIm{1},allMito{1}] = generateMitochondria(numMitoOG,frameSize,micronPerPixel,1,lengthMultiplier,intensityMultiplier);
    close all

    %Generate new mito every frame.
    fusionCheck = 0;
    while fusionCheck == 0
        for frameNum = 2:numFrames

            %Fission, fusion, disappearance, dynamics
            [mitoIm{frameNum},allMito{frameNum}] = generateMitoMovement(allMito{frameNum-1},frameSize,micronPerPixel,secondsPerFrame,frameNum,speedMultiplier,lostChance,fusionChance,fissionChance);

            close all
            %Generate new mito
            numMito = round(normrnd(numMitoOG*percentNewMito,numMitoOG*percentNewMito*0.68));
            if numMito>0
                [newMitoIm{frameNum},mitoSim{frameNum}] = generateMitochondria(numMito,frameSize,0.138,frameNum,1,1);
                close all
                allMito{frameNum} = [allMito{frameNum},mitoSim{frameNum}];
                for newMitoNum = 1:length(mitoSim{frameNum})
                    mitoIm{frameNum}(:,:,end+1) = newMitoIm{frameNum}(:,:,newMitoNum);
                end
            end
        end

        numMito = length(allMito{end});
        track = allMito{1};
        for frameNum = 2:numFrames
            for mitoNum = 1:length(allMito{frameNum})
                if mitoNum > length(allMito{frameNum-1})
                    track(mitoNum) = allMito{frameNum}(mitoNum);
                elseif allMito{frameNum}(mitoNum).fusion
                    track(mitoNum).fusion(end) = allMito{frameNum}(mitoNum).fusion;
                elseif allMito{frameNum}(mitoNum).lost
                    track(mitoNum).lost(end) = allMito{frameNum}(mitoNum).lost;
                elseif allMito{frameNum}(mitoNum).frame == frameNum
                    track(mitoNum) = addToTrack(allMito{frameNum}(mitoNum),track(mitoNum));
                else
                    continue
                end
            end
        end
        if nnz([track.fusion])
            fusionCheck = 1;
        end
    end

    lifetime = zeros(1,numMito);

    for trackNum = 1:numMito
        if find(track(trackNum).fusion)==numFrames
            track(trackNum).fusion(numFrames)=0;
        end
        lifetime(trackNum) = length(track(trackNum).frame);
    end

    dynamics = getDynamics(0.138,5,track);

    comboMitoIm = zeros(frameSize,frameSize,numFrames,'uint8');
    noiseTemp = ones(frameSize,frameSize,'uint8')*noiseMultiplier;
    % noiseTemp = imnoise(noiseTemp,'gaussian');
    for frameNum = 1:numFrames
        comboMitoIm(:,:,frameNum) = sum(mitoIm{frameNum},3);

    %     comboMitoIm(:,:,frameNum) = imgaussfilt(sum(mitoIm{frameNum},3),1);
    %     comboMitoIm(:,:,frameNum) = imnoise(comboMitoIm(:,:,frameNum)+noiseTemp,'Poisson');

    end

    close all
    figure
    for i = 1:size(comboMitoIm,3)
        imagesc(comboMitoIm(:,:,i))
        pause(0.2)
    end

    % User Parameters
    micronPerPixel = 0.138;
    minArea = 0;
    maxArea = 1e4;
    secondsPerFrame = 5;
%     distThreshVals = [0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10,20];
    stdThreshVals = [0,0.1,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,3.5,4];

    numTrials = length(stdThreshVals);

    allAssignments = cell(1,numTrials);
    totalCorrect = zeros(1,numTrials);
    totalAssignments = zeros(1,numTrials);
    percentCorrect = zeros(1,numTrials);
    pValSpeed = zeros(1,numTrials);
    pValLifetime = zeros(1,numTrials);
    timeToc = zeros(1,numTrials);
    aveNNExtrema = zeros(1,numTrials);
    aveNPA = zeros(1,numTrials);

    for trialNum = 1:numTrials
        distThreshMicrons = 1;
        frameThreshSeconds = 15;

        %Preliminary tracking to set weighting
        tic
        weightsPrelim = ones(1,6)*1/6;
        trackPrelim = trackMitochondria(comboMitoIm,comboMitoIm,weightsPrelim,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds,stdThreshVals(trialNum));
        %Get the weights
        [weights,allPerfectTracks] = getTrackingWeights(trackPrelim);

        %Real tracking
        [trackMM,mitoCoM,extra] = trackMitochondria(comboMitoIm,comboMitoIm,weights,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds,stdThreshVals(trialNum));
        %Get the weights
        [weights2,allPerfectTracks2] = getTrackingWeights(trackPrelim);
        timeToc(trialNum) = toc;
        dynamicsMM = getDynamics(0.138,5,trackMM);

        %plot track centroids GT
        close all
        figure
        imagesc(comboMitoIm(:,:,1))
        colormap gray
        caxis([0,255])
        axis image
        hold on
        for i = 1:length(track)
            if length(track(i).frame)>0
                plotWeightedCentroid(track(i))
                plotFissionEvents(track,i)
                plotFusionEvents(track,i)
            end
        end
    
    
        %plot track centroids MM
        % close all
        figure
        imagesc(comboMitoIm(:,:,1))
        colormap gray
        % caxis([-0.05,1]) 
        axis image
        hold on
        for i = 1:length(trackMM)
            if length(trackMM(i).frame)>0
                plotWeightedCentroid(trackMM(i))
                plotFissionEvents(trackMM,i)
                plotFusionEvents(trackMM,i)
            end
        end

        numTracks = length(trackMM);
        numMito = length(track);

        matchTracks = cell(1,numTracks);

        for trackNum = 1:numTracks
            matchTracks{trackNum} = cell(1,length(trackMM(trackNum).frame));
            for frameNum = trackMM(trackNum).frame
                frameIdx = find(trackMM(trackNum).frame==frameNum);
                for mitoNum = 1:numMito
                    GTframeIdx = find(track(mitoNum).frame==frameNum);
                    if isempty(GTframeIdx)
                        continue
                    end
                    if ~iscell(trackMM(trackNum).PixelIdxList)
                        trackMM(trackNum).PixelIdxList = {trackMM(trackNum).PixelIdxList};
                    end
                    if ~iscell(track(mitoNum).PixelIdxList)
                        track(mitoNum).PixelIdxList = {track(mitoNum).PixelIdxList};
                    end
                    if sum(ismember(trackMM(trackNum).PixelIdxList{frameIdx},track(mitoNum).PixelIdxList{GTframeIdx}))
                        matchTracks{trackNum}{frameIdx}(end+1) = mitoNum;
                    end
                end
            end
        end

        lifetimeMM = zeros(1,numTracks);

        correctAssignmentNum = cell(1,numTracks,numTrials);
        for trackNum = 1:numTracks
            numMatches = length(matchTracks{trackNum});
            correctAssignmentNum{trackNum,trialNum} = zeros(3,numMatches-1);
            if numMatches == 1
                correctAssignmentNum{trackNum,trialNum}(1,1) = 0;
                correctAssignmentNum{trackNum,trialNum}(2,1) = trackMM(trackNum).NNExtrema(1);
                correctAssignmentNum{trackNum,trialNum}(3,1) = trackMM(trackNum).NPA(1);
            end
            for matchNum = 2:numMatches
                correctAssignmentNum{trackNum,trialNum}(1,matchNum-1) = sum(ismember(matchTracks{trackNum}{matchNum-1},matchTracks{trackNum}{matchNum}));
                correctAssignmentNum{trackNum,trialNum}(2,matchNum-1) = trackMM(trackNum).NNExtrema(matchNum);
                correctAssignmentNum{trackNum,trialNum}(3,matchNum-1) = trackMM(trackNum).NPA(matchNum);
            end
            lifetimeMM(trackNum) = length(trackMM(trackNum).frame);
            if lifetimeMM(trackNum)<=3
                lifetimeMM(trackNum) = nan;
            end
        end

        allAssignments{trialNum} = [correctAssignmentNum{:,trialNum}];
        totalCorrect(trialNum) = sum(allAssignments{trialNum}(1,:));
        totalAssignments(trialNum) = nnz([track.Area])-length(track);
        percentCorrect(trialNum) = totalCorrect(trialNum)/totalAssignments(trialNum);
        if sum(~isnan([dynamicsMM.speed{:}]))
            pValSpeed(trialNum) = ranksum([dynamicsMM.speed{:}],[dynamics.speed{:}]);
        else
            pValSpeed(trialNum) = 0;
        end
        if sum(~isnan(lifetimeMM))
            pValLifetime(trialNum) = ranksum(lifetimeMM,lifetime);
        else
            pValLifetime(trialNum) = 0;
        end
        aveNNExtrema(trialNum) = nanmean(allAssignments{trialNum}(2,:));
        aveNPA(trialNum) = nanmean(allAssignments{trialNum}(3,:));
        fusionStats(trialNum) = checkFusionStats(track,trackMM,matchTracks);
    end
    comboStats{densityNum} = [percentCorrect;pValSpeed;pValLifetime;aveNNExtrema;aveNPA;timeToc;stdThreshVals];
    fusionStatsCell{densityNum} = fusionStats;
    
end


%%
figure
scatter(allAssignments(3,:),allAssignments(1,:))
%%
%plot perfect track centroids MM
close all
figure
imagesc(comboMitoIm(:,:,1))
colormap gray
axis image
hold on
for i = 1:length(allPerfectTracks2)
    if ~isempty(allPerfectTracks2(i).frame)
        plotWeightedCentroid(allPerfectTracks2(i))
        plotFissionEvents(allPerfectTracks2,i)
        plotFusionEvents(allPerfectTracks2,i)
    end
end

%%
saveTif(comboMitoIm)
%%
simTrack = track;

save('simTrack','simTrack')