function fusionStats = checkFusionStats(track,trackMM,matchTracks)

numTracks = length(trackMM);

TP = 0;
FP = 0;

for trackNum = 1:numTracks
    if ~sum(trackMM(trackNum).fusion)
        continue
    end
    fusionOGIdx = find(trackMM(trackNum).fusion);
    fusionOGFrame = trackMM(trackNum).frame(fusionOGIdx)+1;
    
    fusionWithNum = trackMM(trackNum).fusion(trackMM(trackNum).fusion>0);
    fusionWithIdx = find(trackMM(fusionWithNum).frame>=fusionOGFrame);
%     fusionWithIdx = fusionWithIdxTemp(1);
    
    fusionTrack1 = matchTracks{trackNum}(fusionOGIdx);
    fusionTrack1 = [fusionTrack1{:}];
    fusionTrack2 = matchTracks{fusionWithNum}(fusionWithIdx);
    fusionTrack2 = [fusionTrack2{:}];
    
    check1 = sum(ismember(unique(fusionTrack1),[track(fusionTrack2).fusion]));
    check2 = sum(ismember(unique(fusionTrack2),[track(fusionTrack1).fusion]));
    
    sumCheck = check1+check2;
    
    TP = TP + sumCheck; %top left
    
    if ~sumCheck
        FP = FP + 1; %bottom left
    end
    
end
sumFusionMM = sum(nnz([trackMM.fusion]));

sumFusionGT = sum(nnz([track.fusion])); 
sumNoFusionGT = length(track)-sumFusionGT; 
sumNoFusionMM = length(trackMM)-TP;

FN = max(sumFusionGT-TP,0);
TN = length(track)-FP-FN;

fusionStats.falseNeg = FN;
fusionStats.trueNeg = TN;

fusionStats.truePos = TP;
fusionStats.falsePos = FP;

fusionStats.FPR = FP/(FP+TN);
fusionStats.ACC = (TP+TN)/(TP+TN+FP+FN);

end