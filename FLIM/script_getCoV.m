clear

defaultFileName = fullfile('/Volumes/LEFEBVRE/Mitometer/FLIM', '*.mat');
[fileNames, path] = uigetfile(defaultFileName, 'Select a mat file','MultiSelect','on');

if ~iscell(fileNames)
    fileNames = {fileNames};
end

for fileNum = 1:length(fileNames)
    clear mitoStats
    fullFileName = fullfile(path, fileNames{fileNum});
    load(fullFileName);
    meanfBound(fileNum) = mean(mitoStats.fracB);
    CoVfBound(fileNum) = std(mitoStats.fracB)/mean(mitoStats.fracB);
    
end