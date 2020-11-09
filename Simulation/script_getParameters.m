%% For all control cells
numIms = length(track);

%Analyze only the perfect tracks (tracks with only confident assignments
%and more than 1 frame of analysis)
perfectTrack = cell(1,numIms);
allPerfectTracks = cell(1,numIms);

for imNum = 1:numIms
    numTracks = length(track{imNum});
    for trackNum = 1:numTracks
        if sum(track{imNum}(trackNum).confident)==length(track{imNum}(trackNum).confident) && length(track{imNum}(trackNum).frame)>1
            perfectTrack{imNum}(trackNum) = 1;
        end
    end
    allPerfectTracks{imNum} = track{imNum}(find(perfectTrack{imNum}));
end

%Get all the average areas for every perfect track
allArea = cell(1,numIms);
for imNum = 1:numIms
    numTracks = length(allPerfectTracks{imNum});
    for trackNum = 1:numTracks
        allArea{imNum}(trackNum) = mean(allPerfectTracks{imNum}(trackNum).Area);
    end
end

logArea = log10([allArea{:}]*0.138^2)';
histfit(logArea,20,'normal')
normArea = fitdist(logArea,'normal');
areaMu = normArea.mu;
areaSigma = normArea.sigma;
maxArea = max(logArea);
minArea = min(logArea);
%%
%Get all the average Major Axis Lengths for every perfect track
allMajAx = cell(1,numIms);
for imNum = 1:numIms
    numTracks = length(allPerfectTracks{imNum});
    for trackNum = 1:numTracks
        allMajAx{imNum}(trackNum) = mean(allPerfectTracks{imNum}(trackNum).MajorAxisLength);
    end
end

logMajAx = log10([allMajAx{:}]*0.138)';
histfit(logMajAx,20,'normal')
normMajAx = fitdist(logMajAx,'normal');
majAxMu = normMajAx.mu;
majAxSigma = normMajAx.sigma;
maxMajAx = max(logMajAx);
minMajAx = min(logMajAx);
%%
%Get all the average areas for every perfect track
allMinAx = cell(1,numIms);
for imNum = 1:numIms
    numTracks = length(allPerfectTracks{imNum});
    for trackNum = 1:numTracks
        allMinAx{imNum}(trackNum) = mean(allPerfectTracks{imNum}(trackNum).MinorAxisLength);
    end
end

MinAxRatio = ([allMinAx{:}]./[allMajAx{:}])';
histfit(MinAxRatio,200,'normal')
normMinAxRatio = fitdist(MinAxRatio,'normal');
minAxRatioMu = normMinAxRatio.mu;
minAxRatioSigma = normMinAxRatio.sigma;
maxMinAxRatio = max(MinAxRatio);
minMinAxRatio = min(MinAxRatio);
%%
%Get all the average areas for every perfect track
allInt = cell(1,numIms);
for imNum = 1:numIms
    numTracks = length(allPerfectTracks{imNum});
    for trackNum = 1:numTracks
        allInt{imNum}(trackNum) = mean(allPerfectTracks{imNum}(trackNum).MeanIntensity);
    end
end

logInt = log10([allInt{:}])';
histfit(logInt,20,'normal')
normInt = fitdist(logInt,'normal');
intMu = normInt.mu;
intSigma = normInt.sigma;
%%
allSpeed = cell(1,numIms);
for imNum = 1:numIms
    numTracks = length(allPerfectTracks{imNum});
    for trackNum = 1:numTracks
        allSpeed{imNum}(trackNum) = nanmean(stats{imNum}.speed{trackNum});
    end
end

logSpeed = log10([allSpeed{:}]*0.138/5)';
histfit(logSpeed,20,'normal')
normSpeed = fitdist(logSpeed,'normal');
speedMu = normSpeed.mu;
speedSigma = normSpeed.sigma;
minSpeed = min(logSpeed);
maxSpeed = max(logSpeed);

