clear
addpath(strcat(pwd,'/3D'));

defaultFileName = fullfile('/Volumes/LEFEBVRE/Mitometer/3D', '*.mat');
[fileNames, path] = uigetfile(defaultFileName, 'Select a mat file','MultiSelect','on');

if ~iscell(fileNames)
    fileNames = {fileNames};
end

cellLineName = "Patient2";

for fileNum = 1:length(fileNames)
    fullFileName = fullfile(path, fileNames{fileNum});
    load(fullFileName);
    
    track = trackLengthThreshold(saveTemp.track,saveTemp.optimalThresh);
%     track = saveTemp.allPerfectTracks2;
    dynamics = getDynamics3d(0.1054,10,track,0.45);
    clear mitoStats
    mitoStats = struct();
    for i = 1:length(track)
        mitoStats.MajAx(i) = nanmedian([track(i).MajorAxisLength]);
        mitoStats.MinAx(i) = nanmedian([track(i).MinorAxisLength]);
        mitoStats.ZAx(i) = nanmedian([track(i).ZAxisLength]);
        mitoStats.Sol(i) = nanmedian([track(i).Solidity]);
        mitoStats.Surf(i) = nanmedian([track(i).SurfaceArea]);
        mitoStats.Vol(i) = nanmedian([track(i).Volume]);
        mitoStats.Int(i) = nanmedian([track(i).MeanIntensity]);
        mitoStats.MedianSpeed(i) = nanmedian(dynamics.speed{i});
        mitoStats.MeanSpeed(i) = nanmean(dynamics.speed{i});
        mitoStats.MaxSpeed(i) = max(dynamics.speed{i});
        mitoStats.Directionality(i) = dynamics.displacement{i}(end)/sum(dynamics.distance{i});
        mitoStats.MaxDisplacement(i) = max(dynamics.displacement{i});
        mitoStats.Fission(i) = nnz([track(i).fission]);
        mitoStats.Fusion(i) = nnz([track(i).fusion]);
        
        mitoStats.SampleNum(i) = 12;
        %0 = non-cancerous, 1 = luminal, 2 = basal;
        mitoStats.Invasivity(i) = 2;
        %0 = no, 1 = yes, 2 = unknown;
        mitoStats.ReceptorPos(i) = 0;
        mitoStats.MutTP53(i) = 2;    
        mitoStats.Caucasian(i) = 2; % 0 = Black
        mitoStats.Ductal(i) = 1; % 0 = Adeno, 2 = N/A
        %0 = unknown
        mitoStats.Age(i) = 0;
    end
    
    saveLoc = strcat('/Volumes/LEFEBVRE/Mitometer/3D/','mitoStats_',fileNames{fileNum}(10:end-4),'.mat');
    if ~iscell(saveLoc)
        saveLoc = {saveLoc};
    end
    save(saveLoc{1},'mitoStats');

    aveMito = zeros(1,length(saveTemp.extra.mito));
    for i = 1:length(saveTemp.extra.mito)
        aveMito(i) = length(saveTemp.extra.mito{i});
    end

    aveAll = mean(aveMito);

    cellStats.TotalFission(fileNum) = nnz([track.fission]);
    cellStats.TotalFusion(fileNum) = nnz([track.fusion]);

    cellStats.NumTracks(fileNum) = length(track);

    cellStats.NumMito(fileNum) = nanmean(aveMito);
    cellStats.MeanSpeed(fileNum) = nanmean(mitoStats.MedianSpeed);
    cellStats.MajAx(fileNum) = nanmean(mitoStats.MajAx);
    cellStats.MinAx(fileNum) = nanmean(mitoStats.MinAx);
    cellStats.ZAx(fileNum) = nanmean(mitoStats.ZAx);
    cellStats.Sol(fileNum) = nanmean(mitoStats.Sol);
    cellStats.Surf(fileNum) = nanmean(mitoStats.Surf);
    cellStats.Vol(fileNum) = nanmean(mitoStats.Vol);
    cellStats.Int(fileNum) = nanmean(mitoStats.Int);
    cellStats.MaxSpeed(fileNum) = nanmean(mitoStats.MaxSpeed);
    cellStats.Directionality(fileNum) = nanmean(mitoStats.Directionality);
    cellStats.MaxDisplacement(fileNum) = nanmean(mitoStats.MaxDisplacement);

    cellStats.NumMitoStd(fileNum) = nanstd(aveMito);
    cellStats.MeanSpeedStd(fileNum) = nanstd(mitoStats.MedianSpeed);
    cellStats.MajAxStd(fileNum) = nanstd(mitoStats.MajAx);
    cellStats.MinAxStd(fileNum) = nanstd(mitoStats.MinAx);
    cellStats.ZAxStd(fileNum) = nanstd(mitoStats.ZAx);
    cellStats.SolStd(fileNum) = nanstd(mitoStats.Sol);
    cellStats.SurfStd(fileNum) = nanstd(mitoStats.Surf);
    cellStats.VolStd(fileNum) = nanstd(mitoStats.Vol);
    cellStats.IntStd(fileNum) = nanstd(mitoStats.Int);
    cellStats.MaxSpeedStd(fileNum) = nanstd(mitoStats.MaxSpeed);
    cellStats.DirectionalityStd(fileNum) = nanstd(mitoStats.Directionality);
    cellStats.MaxDisplacementStd(fileNum) = nanstd(mitoStats.MaxDisplacement);

    cellStats.NumMitoCoV(fileNum) = cellStats.NumMitoStd(fileNum)/cellStats.NumMito(fileNum);
    cellStats.MeanSpeedCoV(fileNum) = cellStats.MeanSpeedStd(fileNum)/cellStats.MeanSpeed(fileNum);
    cellStats.MajAxCoV(fileNum) = cellStats.MajAxStd(fileNum)/cellStats.MajAx(fileNum);
    cellStats.MinAxCoV(fileNum) = cellStats.MinAxStd(fileNum)/cellStats.MinAx(fileNum);
    cellStats.ZAxCoV(fileNum) = cellStats.ZAxStd(fileNum)/cellStats.ZAx(fileNum);
    cellStats.SolCoV(fileNum) = cellStats.SolStd(fileNum)/cellStats.Sol(fileNum);
    cellStats.SurfCoV(fileNum) = cellStats.SurfStd(fileNum)/cellStats.Surf(fileNum);
    cellStats.VolCoV(fileNum) = cellStats.VolStd(fileNum)/cellStats.Vol(fileNum);
    cellStats.IntCoV(fileNum) = cellStats.IntStd(fileNum)/cellStats.Int(fileNum);
    cellStats.MaxSpeedCoV(fileNum) = cellStats.MaxSpeedStd(fileNum)/cellStats.MaxSpeed(fileNum);
    cellStats.DirectionalityCoV(fileNum) = cellStats.DirectionalityStd(fileNum)/cellStats.Directionality(fileNum);
    cellStats.MaxDisplacementCoV(fileNum) = cellStats.MaxDisplacementStd(fileNum)/cellStats.MaxDisplacement(fileNum);
end

saveLoc = strcat('/Volumes/LEFEBVRE/Mitometer/3D/','cellStats_',cellLineName,'.mat');
if ~iscell(saveLoc)
    saveLoc = {saveLoc};
end
save(saveLoc{1},'cellStats');

%  for j = 1:length(track)
%     [~,loc] = ismember(track(j).frame,frameArray);
%     out(j,(loc)) = [track(j).(fnames{i})];
% end