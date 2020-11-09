%Images order
% 1 - Intensity
% 2 - Harmonic 1 Phase (degrees)
% 3 - Harmonic 1 Modulation
% 4 - Harmonic 2 Phase (degrees)
% 5 - Harmonic 2 Modulation

function [Images,FileName,PathName]=Open_R64()

[FileName,PathName] = uigetfile('*.R64','MultiSelect','on');
% read .R64 binary file
if ischar(FileName)
    tmp=FileName;
    FileName=cell(1);
    FileName{1}=tmp;
end
L=length(FileName);
Images=cell(1,L);
    for i=1:L
        fileID = fopen([PathName FileName{i}]);
        input = fread(fileID);
        fclose(fileID);
         % decompress file content
         buffer = java.io.ByteArrayOutputStream();
         zlib = java.util.zip.InflaterOutputStream(buffer);
         zlib.write(input, 0, numel(input));
         zlib.close();
         buffer = buffer.toByteArray();
         % read image dimension, number of images, and image data from
         % decompressed buffer
         nimages = 5;
         img=typecast(buffer(5:end), 'single');
         
         sizes=sqrt(length(img)/nimages);
         images = reshape(img, sizes, sizes, nimages);
         Images{i}=double(images);
    end
end
