clear

defaultFileName = fullfile('/Volumes/LEFEBVRE/Mitometer/FLIM', '*.mat');
[fileNames, path] = uigetfile(defaultFileName, 'Select a mat file','MultiSelect','on');

if ~iscell(fileNames)
    fileNames = {fileNames};
end

cellLineName = "BT549";

for fileNum = 1:2:length(fileNames)
    fullFileName = fullfile(path, fileNames{fileNum});
    flimFileName = fullfile(path, fileNames{fileNum+1});
    load(fullFileName);
    load(flimFileName);

    track = trackLengthThreshold(saveTemp.track,20);

    dynamics = getDynamics(0.1758,3.87,track);

    clear mitoStats
    mitoStats = struct();
    for i = 1:length(track)
        mitoStats.fracB(i) = saveFLIM(i).fracB;
        mitoStats.MajAx(i) = nanmedian([track(i).MajorAxisLength]);
        mitoStats.MinAx(i) = nanmedian([track(i).MinorAxisLength]);
%         mitoStats.ZAx(i) = nanmedian([track(i).ZAxisLength]);
        mitoStats.Sol(i) = nanmedian([track(i).Solidity]);
        mitoStats.Peri(i) = nanmedian([track(i).Perimeter]);
        mitoStats.Area(i) = nanmedian([track(i).Area]);
        mitoStats.Int(i) = nanmedian([track(i).MeanIntensity]);
        mitoStats.MedianSpeed(i) = nanmedian(dynamics.speed{i});
        mitoStats.MeanSpeed(i) = nanmean(dynamics.speed{i});
        mitoStats.MaxSpeed(i) = max(dynamics.speed{i});
        mitoStats.Directionality(i) = dynamics.displacement{i}(end)/sum(dynamics.distance{i});
        mitoStats.MaxDisplacement(i) = max(dynamics.displacement{i});
        mitoStats.Fission(i) = nnz([track(i).fission]);
        mitoStats.Fusion(i) = nnz([track(i).fusion]);

        mitoStats.SampleNum(i) = 6;
        %0 = non-cancerous, 1 = luminal, 2 = basal;
        mitoStats.Invasivity(i) = 2;
        %0 = no, 1 = yes, 2 = unknown;
        mitoStats.ReceptorPos(i) = 0;
        mitoStats.MutTP53(i) = 1;    
        mitoStats.Caucasian(i) = 1; % 0 = Black
        mitoStats.Ductal(i) = 1; % 0 = Adeno, 2 = N/A
        %0 = unknown
        mitoStats.Age(i) = 72;
    end

    saveLoc = strcat('/Volumes/LEFEBVRE/Mitometer/FLIM/','mitoStats_',fileNames{fileNum}(10:end-4),'.mat');

    if ~iscell(saveLoc)
        saveLoc = {saveLoc};
    end
    save(saveLoc{1},'mitoStats');

    cellStats.TotalFission(fileNum) = nnz([track.fission]);
    cellStats.TotalFusion(fileNum) = nnz([track.fusion]);

    cellStats.NumTracks(fileNum) = length(track);

    cellStats.fBound(fileNum) = nanmean(mitoStats.fBound);
    cellStats.MeanSpeed(fileNum) = nanmean(mitoStats.MedianSpeed);
    cellStats.MajAx(fileNum) = nanmean(mitoStats.MajAx);
    cellStats.MinAx(fileNum) = nanmean(mitoStats.MinAx);
%     cellStats.ZAx(fileNum) = nanmean(mitoStats.ZAx);
    cellStats.Sol(fileNum) = nanmean(mitoStats.Sol);
    cellStats.Peri(fileNum) = nanmean(mitoStats.Peri);
    cellStats.Area(fileNum) = nanmean(mitoStats.Area);
    cellStats.Int(fileNum) = nanmean(mitoStats.Int);
    cellStats.MaxSpeed(fileNum) = nanmean(mitoStats.MaxSpeed);
    cellStats.Directionality(fileNum) = nanmean(mitoStats.Directionality);
    cellStats.MaxDisplacement(fileNum) = nanmean(mitoStats.MaxDisplacement);

    cellStats.fBoundStd(fileNum) = nanmstd(mitoStats.fBound);
    cellStats.MeanSpeedStd(fileNum) = nanstd(mitoStats.MedianSpeed);
    cellStats.MajAxStd(fileNum) = nanstd(mitoStats.MajAx);
    cellStats.MinAxStd(fileNum) = nanstd(mitoStats.MinAx);
%     cellStats.ZAxStd(fileNum) = nanstd(mitoStats.ZAx);
    cellStats.SolStd(fileNum) = nanstd(mitoStats.Sol);
    cellStats.PeriStd(fileNum) = nanstd(mitoStats.Peri);
    cellStats.AreaStd(fileNum) = nanstd(mitoStats.Area);
    cellStats.IntStd(fileNum) = nanstd(mitoStats.Int);
    cellStats.MaxSpeedStd(fileNum) = nanstd(mitoStats.MaxSpeed);
    cellStats.DirectionalityStd(fileNum) = nanstd(mitoStats.Directionality);
    cellStats.MaxDisplacementStd(fileNum) = nanstd(mitoStats.MaxDisplacement);

    cellStats.fBound(fileNum) = cellStats.fBoundStd(fileNum)/cellStats.fBound(fileNum);
    cellStats.MeanSpeedCoV(fileNum) = cellStats.MeanSpeedStd(fileNum)/cellStats.MeanSpeed(fileNum);
    cellStats.MajAxCoV(fileNum) = cellStats.MajAxStd(fileNum)/cellStats.MajAx(fileNum);
    cellStats.MinAxCoV(fileNum) = cellStats.MinAxStd(fileNum)/cellStats.MinAx(fileNum);
%     cellStats.ZAxCoV(fileNum) = cellStats.ZAxStd(fileNum)/cellStats.ZAx(fileNum);
    cellStats.SolCoV(fileNum) = cellStats.SolStd(fileNum)/cellStats.Sol(fileNum);
    cellStats.PeriCoV(fileNum) = cellStats.PeriStd(fileNum)/cellStats.Peri(fileNum);
    cellStats.AreaCoV(fileNum) = cellStats.AreaStd(fileNum)/cellStats.Area(fileNum);
    cellStats.IntCoV(fileNum) = cellStats.IntStd(fileNum)/cellStats.Int(fileNum);
    cellStats.MaxSpeedCoV(fileNum) = cellStats.MaxSpeedStd(fileNum)/cellStats.MaxSpeed(fileNum);
    cellStats.DirectionalityCoV(fileNum) = cellStats.DirectionalityStd(fileNum)/cellStats.Directionality(fileNum);
    cellStats.MaxDisplacementCoV(fileNum) = cellStats.MaxDisplacementStd(fileNum)/cellStats.MaxDisplacement(fileNum);
end

saveLoc = strcat('/Volumes/LEFEBVRE/Mitometer/FLIM/','cellStats_',cellLineName,'.mat');

if ~iscell(saveLoc)
    saveLoc = {saveLoc};
end
save(saveLoc{1},'cellStats');

%  for j = 1:length(track)
%     [~,loc] = ismember(track(j).frame,frameArray);
%     out(j,(loc)) = [track(j).(fnames{i})];
% end