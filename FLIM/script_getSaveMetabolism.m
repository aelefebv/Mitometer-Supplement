% clear
% close all
% clc
addpath(strcat(pwd,'/FLIM'))
sprintf("Choose the NADH R64 files in order.")
NADHin = Save_R64();

strcat(NADHin(end).path,NADHin(end).name)

% sprintf("Choose the corresponding .mat file.")
% 
% defaultFileName = fullfile('/Volumes/LEFEBVRE/Mitometer/FLIM', '*.mat');
% [fileName, path] = uigetfile(defaultFileName, 'Select a mat file','MultiSelect','off');
% 
% if ~iscell(fileName)
%     fileName = {fileName};
% end
% 
% fullFileName = fullfile(path, fileName{1});
% load(fullFileName);
%%
addpath(strcat(pwd,'/FLIM'))

for i = 1:length(NADHin)
    NADHint(:,:,i) = NADHin(i).intensity;
    NADHg(:,:,i) = NADHin(i).g;
    NADHs(:,:,i) = NADHin(i).s;
end
% NADH = im2uint8(NADHint);

clear FLIMtrack allG allS

kept = trackLengthThreshold(saveTemp.track,5);

for i = 1:length(kept)
    for j = 1:length(kept(i).frame)
        tempInt = NADHint(:,:,kept(i).frame(j));
        tempG = NADHg(:,:,kept(i).frame(j));
        tempS = NADHs(:,:,kept(i).frame(j));
        pixelList = kept(i).PixelIdxList{j};
        FLIMtrack(i).int{j} = tempInt(pixelList);
        FLIMtrack(i).g{j} = tempG(pixelList(FLIMtrack(i).int{j}>0));
        FLIMtrack(i).s{j} = tempS(pixelList(FLIMtrack(i).int{j}>0));
%         FLIMtrack(i).allG = [];
%         FLIMtrack(i).allS = [];
%         FLIMtrack(i).meanG = [];
%         FLIMtrack(i).meanS = [];
    end
end

for i = 1:length(FLIMtrack)
    allG{i} = [];
    allS{i} = [];
    for j = 1:length(FLIMtrack(i).g)
        allG{i} = [allG{i}; FLIMtrack(i).g{j}];
        allS{i} = [allS{i}; FLIMtrack(i).s{j}];
    end
    FLIMtrack(i).allG = allG{i};
    FLIMtrack(i).allS = allS{i};

    meanG(i) = mean(allG{i});
    meanS(i) = mean(allS{i});
    
    FLIMtrack(i).meanG = meanG(i);
    FLIMtrack(i).meanS = meanS(i);
    
    [frac1,frac2,frac3] = fractionCalc(meanS(i),meanG(i),0.4E-9,3.2E-9,80000000);
    
    fracF = frac1/(frac1+frac2);
    fracB = frac2/(frac1+frac2);
    
    FLIMtrack(i).fracF = fracF;
    FLIMtrack(i).fracB = fracB;
    
end
%% label mitochondria with their fractionBound

labelIm=zeros(253,253);
for i = 1:length(kept)
tempLabel = plotLables(kept,extra.L,i)*FLIMtrack(i).fracB;
labelIm = labelIm + tempLabel;
end
imagesc(labelIm(:,:,1))
colormap gray
caxis([0,1])
axis image off
figure
imagesc(ImOG(:,:,1))
axis image off
colormap gray
% saveFLIM = FLIMtrack;  
% saveLoc = strcat(path,fileName{1}(1:end-4),'_FLIM.mat');
% if ~iscell(saveLoc)
%     saveLoc = {saveLoc};
% end
% save(saveLoc{1},'saveFLIM');
