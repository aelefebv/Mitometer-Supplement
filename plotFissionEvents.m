function plotFissionEvents(track,trackNum)
%Enter a specific track (i.e. track(#)) whose weighted centroids you wish
%to see.

if nnz(track(trackNum).fission)>0
    fissionTrack = track(trackNum).fission(track(trackNum).fission>0);
    fissionFrame = track(trackNum).frame(track(trackNum).fission>0)-1;
    OGTrackFrame = track(fissionTrack).frame(track(fissionTrack).frame==fissionFrame);
    while isempty(OGTrackFrame)
        fissionFrame = fissionFrame-1;
        OGTrackFrame = track(fissionTrack).frame(track(fissionTrack).frame==fissionFrame);
    end

    OGTrackCentroid1 = track(fissionTrack).WeightedCentroid(find(track(fissionTrack).frame == OGTrackFrame)*2-1);
    OGTrackCentroid2 = track(fissionTrack).WeightedCentroid(find(track(fissionTrack).frame == OGTrackFrame)*2);
    connect1 = [OGTrackCentroid1,track(trackNum).WeightedCentroid(1)];
    connect2 = [OGTrackCentroid2,track(trackNum).WeightedCentroid(2)];

%     plot(connect1,connect2,':x','LineWidth',2,'Color','#7E2F8E')
    plot(connect1,connect2,'m.','MarkerSize',20)
end


end