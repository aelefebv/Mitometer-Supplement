function [FlimIntensity,FlimG,FlimS] = SaveIntGS(FlimIm,saveFiles)
%Input your FLIM structure. Put 1 as saveFiles if you want to save the
%intensity, G, and S tif files in the original folder.

numFrames = length(FlimIm);

FlimIntensity = zeros(size(FlimIm(1).intensity,1),size(FlimIm(1).intensity,2),numFrames);
FlimG = zeros(size(FlimIm(1).intensity,1),size(FlimIm(1).intensity,2),numFrames);
FlimS = zeros(size(FlimIm(1).intensity,1),size(FlimIm(1).intensity,2),numFrames);

for i = 1:numFrames
    FlimIntensity(:,:,i) = FlimIm(i).intensity;
    FlimG(:,:,i) = FlimIm(i).g;
    FlimS(:,:,i) = FlimIm(i).s;
end

if saveFiles
    saveTif(im2uint8(FlimIntensity),FlimIm(1).path,strcat('int_',FlimIm(1).name,'.tif'));
    saveTif(FlimG,FlimIm(1).path,strcat('g_',FlimIm(1).name,'.tif'));
    saveTif(FlimS,FlimIm(1).path,strcat('s_',FlimIm(1).name,'.tif'));
end

end