function [mitoIm,mitoSim] = generateSpecificMitochondria(majorLength,extremaDist,vert,noiseLevel,imageSize)

if nargin<5
    imageSize = 512;
end
micronPerPixel = 0.138;
frameNum = 1;
numMito = 2;

im = ones(imageSize,imageSize);

mainFig = figure;
set(mainFig, 'MenuBar', 'none')
set(mainFig, 'Toolbar', 'none')

mitoIm = zeros(imageSize,imageSize,numMito,'uint8');
for mitoNum = 1:numMito
    
    imshow(im)
    axis image off
    hold on

    %Override majorAxis
    a = majorLength/micronPerPixel/2;

    %Override minorAxis
    b = 1.1555/micronPerPixel/4;

    %random initial position (constrained)
    x0 = imageSize/2;
    
    if vert == 1
        rotate = pi/2;
        lengthCheck = a;
    elseif vert == 2
        rotate = pi/3;
        lengthCheck = a;
    else
        rotate = 0;
        lengthCheck = b;
    end
    
    if mitoNum == 1
        y0 = imageSize/2 + lengthCheck + extremaDist/2;
    else
        y0 = imageSize/2 - lengthCheck - extremaDist/2;
    end

    %rotate mitochondrion

    
    %fill in the mitochondria and plot them.
    drawMito([x0,y0],a*2,b*2,rotate)

    %Make the plot into an image.
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    rect = [pos(1),pos(2),pos(3),pos(4)];
    f = getframe(gcf,rect);
    [X,~] = frame2im(f);

    %inverse colors for tracking.
    tempMitoIm = 1-X(:,:,1);

    %override Intensity
    intensity = max(1,255-noiseLevel);
    
    mitoIm(:,:,mitoNum) = tempMitoIm*intensity;

    %get mito properties (will be slightly different than generated and 
    %plotted stats)
    CCSim = bwconncomp(mitoIm(:,:,mitoNum));
    mito = regionprops(CCSim,mitoIm(:,:,mitoNum),'Area','Centroid','Eccentricity','MaxFeretProperties','Extrema','MeanIntensity','WeightedCentroid','MajorAxisLength','MinorAxisLength','Perimeter','PixelIdxList');

    
    %Correct the angles of the mitochondria and assign NN, frame, and
    %initiliaze other parameters.
    mito.MaxFeretAngle = angle0to180(180-mito.MaxFeretAngle);
    mito.NN =0;
    mito.NNExtrema = 0;
    mito.frame = frameNum;
    mito.NPA = 0;
    mito.confident = 1;
    mito.lost = 0;
    mito.fission = 0;
    mito.fusion = 0;
    mito.a = a;
    mito.b = b;
    mito.rotate = rotate;
    mito.intensity = intensity;
    mito.trueCoords = [x0,y0];
    
    if numMito == 1
        mitoSim = mito;
    else
        mitoSim(mitoNum) = mito;
    end
    
    hold off
end
close
end
