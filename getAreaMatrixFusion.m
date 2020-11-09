function areaMatrix = getAreaMatrixFusion(lostTrack,track,frameNum)

numTracks = length(track);

%initialize
areaMatrix = zeros(1,numTracks);

for trackNum = 1:numTracks
    frameIdx = find(track(trackNum).frame==frameNum);
    if frameIdx
        areaMatrix(1,trackNum) = -(lostTrack.Area(end)-track(trackNum).Area(frameIdx))./(track(trackNum).Area(frameIdx));
    else
        areaMatrix(1,trackNum) = Inf;
    end
end
end
