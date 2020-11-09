function extremaMatrix = getExtremaMatrix(mito,track)

numMito = length(mito);
numTracks = length(track);

%initialize
extremaMatrix = zeros(numMito,numTracks);

for mitoNum = 1:numMito
    for trackNum = 1:numTracks
        extremaMatrix(mitoNum,trackNum) = min(sqrt( (mito(mitoNum).Extrema(:,1)-track(trackNum).Extrema(:,end-1)).^2 + (mito(mitoNum).Extrema(:,2)-track(trackNum).Extrema(:,end)).^2 ));
    end
end

end
