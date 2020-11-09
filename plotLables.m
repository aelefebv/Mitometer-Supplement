function labelImage = plotLables(track,L,mitoNum)
%Enter a specific track (i.e. track(#)) whose weighted centroids you wish
%to see.

numFrames = max([track.frame]);
labelImage = zeros(size(L{1},1),size(L{1},2),numFrames);

for frameNum = 1:numFrames
    idx = find(track(mitoNum).frame==frameNum);
    if idx
        labelImage(:,:,frameNum) = (L{frameNum} == track(mitoNum).label(idx(end)));
    end
end

end