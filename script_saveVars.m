groupName = 'high';
holdIms = high;
micronPerPixel = 0.138;
secondsPerFrame = 5;
thresholdNum = 3;

% keptTracks = trackLengthThreshold(holdIms(1).track,holdIms(1).optimalThresh);
keptTracks = trackLengthThreshold(holdIms(1).track,thresholdNum);
for i = 2:length(holdIms)
    holdTracks = trackLengthThreshold(holdIms(i).track,thresholdNum);
    keptTracks = [keptTracks,holdTracks];
end

dynamics = getDynamics(micronPerPixel,secondsPerFrame,keptTracks);
for i = 1:length(keptTracks)
keptTracks(i).speed = dynamics.speed{i};
keptTracks(i).distance = dynamics.distance{i};
keptTracks(i).displacement = dynamics.displacement{i};
end

maxFrames = max([keptTracks.frame]);
frameArray = 1:maxFrames;
numTracks = length(keptTracks);
fnames = {'speed','distance','displacement'};


savePath = uigetdir('Select a folder to save your files. Text files will be named "date and time"+"_group name"+"_parameter name"');
for i = 1:length(fnames)
    out = nan(numTracks,maxFrames);
    for j = 1:length(keptTracks)
        [~,loc] = ismember(keptTracks(j).frame,frameArray);
        out(j,(loc)) = [keptTracks(j).(fnames{i})];
    end
    writematrix(out,strcat(savePath,"/",datestr(now,'yyyymmddHHMMSS_'),groupName,'_',fnames{i},'.txt'))
end

