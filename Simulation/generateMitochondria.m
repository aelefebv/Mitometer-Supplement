function [mitoIm,mitoSim] = generateMitochondria(numMito,imageSize,micronPerPixel,frameNum,lengthMultiplier,intensityMultiplier)

im = ones(imageSize,imageSize);

mainFig = figure;
set(mainFig, 'MenuBar', 'none')
set(mainFig, 'Toolbar', 'none')

mitoIm = zeros(imageSize,imageSize,numMito,'uint8');
for mitoNum = 1:numMito
    
    imshow(im)
    axis image off
    hold on
    
    %generate Major Axis Length from lognormal distribution
    mu = 0.518344;
    sigma = 0.336479;
    MajAxlog = normrnd(mu,sigma);
    while MajAxlog > 1.4885 || MajAxlog < -0.0990
        MajAxlog = normrnd(mu,sigma);
    end

    %Major axis from the log normal distribution above
    a = 10^MajAxlog/micronPerPixel/4*lengthMultiplier;

    % %generate Minor Axis Length Ratio from lognormal distribution
    mu = 0.446964;
    sigma = 0.18121;
    MinAxRatio = normrnd(mu,sigma);
    while MinAxRatio > 0.9808 || MinAxRatio < 0.0961
        MinAxRatio = normrnd(mu,sigma);
    end

    %Minor axis is the median of control values
    
    %Override minorAxis
    b = 1.1555/micronPerPixel/4;

    %random initial position (constrained)
    x0 = unifrnd(imageSize/3,imageSize-imageSize/3);
    y0 = unifrnd(imageSize/3,imageSize-imageSize/3);

    %rotate mitochondrion
    rotate = -unifrnd(0,pi)-(pi/2);

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

    %Generate intensity from a log normal (from control data)
    mu = 1.20268;
    sigma = 0.254312;
    intensityLog = normrnd(mu,sigma);
    intensity = 10^intensityLog*intensityMultiplier;
    while intensity > 255 || intensity < 1
        intensityLog = normrnd(mu,sigma);
        intensity = 10^intensityLog*intensityMultiplier;
    end

    %override Intensity
    intensity = normrnd(100,1);
    
%     mitoIm(:,:,mitoNum) = imgaussfilt(tempMitoIm*intensity,0.25/micronPerPixel);
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
    mito.overlap1 = 0;
    
    if numMito == 1
        mitoSim = mito;
    else
        mitoSim(mitoNum) = mito;
    end
    
    hold off
end
close
end
