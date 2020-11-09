numTracks = length(track);

for trackNum = 1:numTracks
    numFrames = length(track(trackNum).frame);
    for frameNum = 1:numFrames
        frameIdx = track(trackNum).frame(frameNum);
        pixelIDs = track(trackNum).PixelIdxList{frameNum};
        
        trackPhase = NADH(frameIdx).phase(pixelIDs);
        trackMod = NADH(frameIdx).mod(pixelIDs);
        trackInt = NADH(frameIdx).intensity(pixelIDs);
        trackInt(trackInt<=2) = 0;
        trackSt(trackNum).phaseSum(frameNum) = sum(trackPhase.*trackInt);
        trackSt(trackNum).modSum(frameNum) = sum(trackMod.*trackInt);
        trackSt(trackNum).intSum(frameNum) = sum(trackInt);

    end
    trackSt(trackNum).phase = sum(trackSt(trackNum).phaseSum)/sum(trackSt(trackNum).intSum);
    trackSt(trackNum).mod = sum(trackSt(trackNum).modSum)/sum(trackSt(trackNum).intSum);
    
    trackSt(trackNum).g = trackSt(trackNum).mod.*cosd(trackSt(trackNum).phase);
    trackSt(trackNum).s = trackSt(trackNum).mod.*sind(trackSt(trackNum).phase);
    [f1,f2,f3] = fractionCalc(trackSt(trackNum).s,trackSt(trackNum).g,0.4e-9,3.2e-9,80000000);
    
    trackSt(trackNum).fBound = f2/(f1+f2);
end