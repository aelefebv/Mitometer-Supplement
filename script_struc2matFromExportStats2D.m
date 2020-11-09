clear

defaultFileName = fullfile('/Volumes/LEFEBVRE/Mitometer/2D', '*.mat');
[fileNames, path] = uigetfile(defaultFileName, 'Select a mat file','MultiSelect','on');

if ~iscell(fileNames)
    fileNames = {fileNames};
end
allOut = [];
for fileNum = 1:length(fileNames)
    fullFileName = fullfile(path, fileNames{fileNum});
    load(fullFileName);
    allOut = [allOut,cell2mat(struct2cell(mitoStats))];
end
allOutHeaders = {'MajAx', 'MinAx', 'Sol', 'Peri', 'Area', 'Int', 'MedianSpeed', 'MeanSpeed', 'MaxSpeed', 'Directionality', 'MaxDisplacement', 'Fission', 'Fusion', 'SampleNum', 'Invasivity', 'ReceptorPos', 'MutTP53', 'Caucasian', 'Ductal', 'Age'};
writeOut = [allOutHeaders;num2cell(allOut')];
writecell(writeOut,strcat(path,"/",'allMito','.txt'));
