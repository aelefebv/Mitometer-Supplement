function saveFLIMtif(name)

flimIm = Save_R64();
imInt = zeros(size(flimIm(1).intensity,1),size(flimIm(1).intensity,2),length(flimIm));
for i = 1:length(flimIm)
imInt(:,:,i) = flimIm(i).intensity;
end
Im = uint8(imInt);

if ~iscell(name)
    name = {name};
end

outputFileName = strcat('/Volumes/LEFEBVRE/Mitometer/FLIM/',"combined_",name{1},"_",flimIm(1).name,".tif");
for iii=1:size(Im,3)
    for iv = 1:size(Im,4)
        imwrite(Im(:,:,iii,iv), outputFileName, 'WriteMode', 'append',  'Compression','none');
    end
end

end