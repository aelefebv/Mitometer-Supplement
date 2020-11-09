function saveTif(Im,path,inputName)
if nargin==1
    [inputName,path,indx] = uiputfile('*.tif','Choose a folder in which to save selected image(s)');
    if indx == 0
        sprintf('User canceled.')
        return
    end
end
outputFileName = strcat(path,inputName);
if exist(outputFileName,"file")
    overwritePrompt = sprintf('File already exists. Overwrite?');
    overwriteButton = questdlg(overwritePrompt, 'Overwrite?', 'Overwrite', 'Cancel', 'Overwrite');

    if strcmpi(overwriteButton, 'Overwrite')
        delete(outputFileName)
        for iii=1:size(Im,3)
            imwrite(Im(:, :, iii), outputFileName, 'WriteMode', 'append',  'Compression','none');
        end
    end
else
    for iii=1:size(Im,3)
            imwrite(Im(:, :, iii), outputFileName, 'WriteMode', 'append',  'Compression','none');
    end
end
end
