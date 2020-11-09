function fissionStats = checkFissionStats(track,trackMM,matchTracks)

numTracks = length(trackMM);

TP = 0;
FP = 0;

for trackNum = 1:numTracks
    if ~sum(trackMM(trackNum).fission)
        continue
    end
    fissionOGIdx = find(trackMM(trackNum).fission);
    fissionOGFrame = trackMM(trackNum).frame(fissionOGIdx)-1;
    
    fissionWithNum = trackMM(trackNum).fission(trackMM(trackNum).fission>0);
    fissionWithIdx = find(trackMM(fissionWithNum).frame<=fissionOGFrame);
%     fissionWithIdx = fissionWithIdxTemp(1);
    
    fissionTrack1 = matchTracks{trackNum}(fissionOGIdx);
    fissionTrack1 = [fissionTrack1{:}];
    fissionTrack2 = matchTracks{fissionWithNum}(fissionWithIdx);
    fissionTrack2 = [fissionTrack2{:}];
    
    check1 = sum(ismember(unique(fissionTrack1),[track(fissionTrack2).fission]));
    check2 = sum(ismember(unique(fissionTrack2),[track(fissionTrack1).fission]));
    
    sumCheck = check1+check2;
    
    TP = TP + sumCheck; %top left
    
    if ~sumCheck
        FP = FP + 1; %bottom left
    end
    
end
sumFissionMM = sum(nnz([trackMM.fission]));

sumFissionGT = sum(nnz([track.fission])); 
sumNoFissionGT = length(track)-sumFissionGT; 
sumNoFissionMM = length(trackMM)-TP;

FN = max(sumFissionGT-TP,0);
TN = length(track)-FP-FN;

fissionStats.falseNeg = FN;
fissionStats.trueNeg = TN;

fissionStats.truePos = TP;
fissionStats.falsePos = FP;


fissionStats.FPR = FP/(FP+TN);
fissionStats.ACC = (TP+TN)/(TP+TN+FP+FN);



end