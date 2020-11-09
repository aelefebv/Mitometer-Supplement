for i = 1:length(searchTimeSave)
    im = searchTimeSave(i);
    for j = 1:length(im.searchTimeVal)
        
        trackVal = im.searchTimeVal(j);
%         [~,perfTrack,~] = getTrackingWeights(trackVal.track);
        dyn = getDynamics(1,1,trackVal.track);
        speedStatsAll = getDynamicStats(dyn.speed);
        speedStats = speedStatsAll.mean;
        pvalSpeed(1:3,1) = nan;
        if j>1
            pvalSpeed(i,j) = ranksum(speedStats,speedStatsPre);
        end
        speedStatsPre = speedStats;
%         meanMeanSpeed(i,j) = nanmean(speedStats.mean);
%         trackNum(i,j) = length(trackVal.track);
%         perfTrackNum(i,j) = length(perfTrack);

%         for k = 1:length(trackVal.track)
%             trk = trackVal.track(k);
%             CoVarea(k) = nanstd(trk.Area)/nanmean(trk.Area);
%             CoVint(k) = nanstd(trk.MeanIntensity)/nanmean(trk.MeanIntensity);
%             CoVsol(k) = nanstd(trk.Solidity)/nanmean(trk.Solidity);
%         end
%         meanCoVarea(i,j) = nanmean(CoVarea(CoVarea~=0));
%         meanCoVint(i,j) = nanmean(CoVint(CoVint~=0));
%         meanCoVsol(i,j) = nanmean(CoVsol(CoVsol~=0));
    end
end