function plotExtrema(track)
%Enter a specific track (i.e. track(#)) whose extremas you wish to see.
all1 = track.Extrema(:,1:2:end);
all2 = track.Extrema(:,2:2:end);
plot(all1(:),all2(:))
end