function [newMitoIm,newMito] = generateMitoMovement(mito,frameSize,micronPerPixel,secondsPerFrame,frameNum,speedMultiplier,lostChance,fusionChance,fissionChance)

im = ones(frameSize,frameSize);

mainFig = figure;
set(mainFig, 'MenuBar', 'none')
set(mainFig, 'Toolbar', 'none')

newMitoIm = zeros(frameSize,frameSize,length(mito),'uint8');

newMito = mito;

travelDistanceAll = zeros(1,length(mito));

randMitoList = (1:length(mito));
randMitoList = randMitoList(randperm(length(randMitoList)));

for mitoNum = randMitoList
% while mitoNum < length(mito)
%     mitoNum = mitoNum + 1;
    %skip if already disappeared or fused
    imshow(im)
    axis image off
    hold on
    
    if mito(mitoNum).lost(end) || mito(mitoNum).fusion(end)
        
        continue
    end

    %Chance to disappear
    if rand<lostChance
        mito(mitoNum).lost(end) = 1;
        newMito(mitoNum).lost(end) = 1;
        continue
    end
    
    %generate speed from log normal distribution
    mu = -1.45484;
    sigma = 0.367769;
    travelDistanceLog = normrnd(mu,sigma);
    while (travelDistanceLog < -2.9010) || (travelDistanceLog > -0.4621)
        travelDistanceLog = normrnd(mu,sigma);
    end
    
    travelDistance = 10^travelDistanceLog*secondsPerFrame/micronPerPixel*speedMultiplier;
    %chance for backwards movement
    if rand < 0.2239
        travelDistance = -travelDistance;
    end
    travelDistanceAll(mitoNum) = travelDistance;
    
    %rotate mitochondrion
    rotate = normrnd(mito(mitoNum).rotate(end),deg2rad(1)*secondsPerFrame);

    %move mitochondrion in direction of rotation
    newX = mito(mitoNum).trueCoords(1)+travelDistance*cos(rotate);
    newY = mito(mitoNum).trueCoords(2)+travelDistance*sin(rotate);
    
    %fill in the mitochondria and plot them.
    drawMito([newX,newY],normrnd(mito(mitoNum).a*2,0.1),normrnd(mito(mitoNum).b*2,0.1),rotate)

    
    %Make the plot into an image.
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    rect = [pos(1),pos(2),pos(3),pos(4)];
    f = getframe(gcf,rect);
    [X,~] = frame2im(f);
    
    %inverse colors for tracking.
    tempMitoIm = 1-X(:,:,1);
%     newMitoIm(:,:,mitoNum) = imgaussfilt(tempMitoIm*mito(mitoNum).intensity,0.25/micronPerPixel);
    newMitoIm(:,:,mitoNum) = tempMitoIm*normrnd(mito(mitoNum).intensity,0.1);
    %     
%     close all
%     imagesc(mitoIm(:,:,1)+newMitoIm(:,:,1))
%     hold on
%     plot([mito(1).WeightedCentroid(1),newX],[mito(1).WeightedCentroid(2),newY])
    
    %get mito properties (will be slightly different than generated and 
    %plotted stats)
    CCNew = bwconncomp(newMitoIm(:,:,mitoNum));
    newMitoTemp = regionprops(CCNew,newMitoIm(:,:,mitoNum),'Area','Centroid','Eccentricity','MaxFeretProperties','Extrema','MeanIntensity','WeightedCentroid','MajorAxisLength','MinorAxisLength','Perimeter','PixelIdxList');
    
    %Correct the angles of the mitochondria and assign NN, frame, and
    %initiliaze other parameters.
    newMitoTemp.MaxFeretAngle = angle0to180(180-newMitoTemp.MaxFeretAngle);
    newMitoTemp.NN =0;
    newMitoTemp.NNExtrema = 0;
    newMitoTemp.frame = frameNum;
    newMitoTemp.NPA = 0;
    newMitoTemp.confident = 1;
    newMitoTemp.lost = 0;
    newMitoTemp.fission = 0;
    newMitoTemp.fusion = mito(mitoNum).fusion;
    newMitoTemp.a = mito(mitoNum).a;
    newMitoTemp.b = mito(mitoNum).b;
    newMitoTemp.rotate = rotate;
    newMitoTemp.intensity = mito(mitoNum).intensity;
    newMitoTemp.trueCoords = [newX,newY];
    newMitoTemp.overlap1 = mito(mitoNum).overlap1(end);
   
    newMito(mitoNum) = newMitoTemp;
    
    hold off
end

overlap = zeros(length(newMito),length(newMito));

%check if there is overlap
for mitoNum = 1:length(newMito)
    if mito(mitoNum).lost(end) || mito(mitoNum).fusion(end)
        continue
    end
    for overlapCheck = mitoNum:length(newMito)
        if overlapCheck == mitoNum
            continue
        end
        overlap(mitoNum,overlapCheck) = sum(ismember(mito(mitoNum).PixelIdxList,mito(overlapCheck).PixelIdxList));
    end
end
%chance for fusion if there is overlap

for mitoNum = randMitoList
    if mito(mitoNum).lost(end) || mito(mitoNum).fusion(end)
        continue
    end
    if ~nnz(overlap(mitoNum,:))
        continue
    end
    for overlapNum = find(overlap(mitoNum,:))
        if mito(mitoNum).overlap1 == 1 || ~overlap(mitoNum,overlapNum) || sum(mito(mitoNum).fusion) || sum(mito(overlapNum).fusion) || sum(mito(overlapNum).fusion)
            continue
        end
        if mito(overlapNum).lost(end)
            continue
        end
        if mito(mitoNum).frame(end) == 1
            mito(overlapNum).overlap1(:) = 1;
            newMito(overlapNum).overlap1 = 1;
            mito(mitoNum).overlap1(:) = 1;
            newMito(mitoNum).overlap1 = 1;
            continue
        end
        temp = rand;
%         if temp>=fusionChance
%             newMito(mitoNum) = mito(mitoNum);
%             newMito(mitoNum).lost(end) = 1;
%             newMitoIm(:,:,mitoNum) = zeros(frameSize,frameSize);
%         elseif temp<fusionChance && ~newMito(mitoNum).fusion && (newMito(mitoNum).a+newMito(overlapNum).a)<111.5812 %max mito length
        if temp<fusionChance && ~newMito(mitoNum).fusion(end) && ~newMito(overlapNum).fusion(end) && ~mito(overlapNum).fusion(end) && ~mito(mitoNum).fusion(end) && (newMito(mitoNum).a+newMito(overlapNum).a)<111.5812 %max mito length

            imshow(im)
            axis image off
            hold on
            
            newMito(mitoNum) = mito(mitoNum);
            newMito(mitoNum).fusion(end) = overlapNum;
            overlap(mitoNum,:) = 0;
            overlap(:,overlapNum) = 0;
            newMitoIm(:,:,mitoNum) = zeros(frameSize,frameSize);
            
            newMito(overlapNum) = simulateFusion(newMito(mitoNum),newMito(overlapNum));
            
            newX = newMito(overlapNum).trueCoords(1);
            newY = newMito(overlapNum).trueCoords(2);
            %fill in the mitochondria and plot them.
            drawMito([newX,newY],normrnd(newMito(overlapNum).a*2,0.1),normrnd(newMito(overlapNum).b*2,0.1),newMito(overlapNum).rotate)

            %Make the plot into an image.
            ax = gca;
            ax.Units = 'pixels';
            pos = ax.Position;
            rect = [pos(1),pos(2),pos(3),pos(4)];
            f = getframe(gcf,rect);
            [X,~] = frame2im(f);

            %inverse colors for tracking.
            tempMitoIm = 1-X(:,:,1);
%             newMitoIm(:,:,overlapNum) = imgaussfilt(tempMitoIm*newMito(overlapNum).intensity,0.25/micronPerPixel);
            newMitoIm(:,:,overlapNum) = tempMitoIm*normrnd(newMito(overlapNum).intensity,0.1);
            %     
        %     close all
        %     imagesc(mitoIm(:,:,1)+newMitoIm(:,:,1))
        %     hold on
        %     plot([mito(1).WeightedCentroid(1),newX],[mito(1).WeightedCentroid(2),newY])

            %get mito properties (will be slightly different than generated and 
            %plotted stats)
            CCNew = bwconncomp(newMitoIm(:,:,overlapNum));
            newMitoTemp = regionprops(CCNew,newMitoIm(:,:,overlapNum),'Area','Centroid','Eccentricity','MaxFeretProperties','Extrema','MeanIntensity','WeightedCentroid','MajorAxisLength','MinorAxisLength','Perimeter','PixelIdxList');

            %Correct the angles of the mitochondria and assign NN, frame, and
            %initiliaze other parameters.
            newMitoTemp.MaxFeretAngle = angle0to180(180-newMitoTemp.MaxFeretAngle);
            newMitoTemp.NN = newMito(overlapNum).NN;
            newMitoTemp.NNExtrema = newMito(overlapNum).NNExtrema;
            newMitoTemp.frame = newMito(overlapNum).frame;
            newMitoTemp.NPA = newMito(overlapNum).NPA;
            newMitoTemp.confident = 1;
            newMitoTemp.lost = newMito(overlapNum).lost;
            newMitoTemp.fission = newMito(overlapNum).fission;
            newMitoTemp.fusion = newMito(overlapNum).fusion;
            newMitoTemp.a = newMito(overlapNum).a;
            newMitoTemp.b = newMito(overlapNum).b;
            newMitoTemp.rotate = newMito(overlapNum).rotate;
            newMitoTemp.intensity = newMito(overlapNum).intensity;
            newMitoTemp.trueCoords = [newX,newY];
            newMitoTemp.overlap1 = mito(mitoNum).overlap1;


            newMito(overlapNum) = newMitoTemp;
            
            hold off
            continue
        end
    end
end

%chance to undergo fission
for mitoNum = randMitoList
    if mito(mitoNum).lost(end) || mito(mitoNum).fusion(end)
        continue
    end
    if frameNum < 3
        continue
    end
    if rand>=fissionChance || newMito(mitoNum).fission || newMito(mitoNum).lost || newMito(mitoNum).a/2<newMito(mitoNum).b
        continue
    end
    
    imshow(im)
    axis image off
    hold on
            
    fissionMito = newMito(mitoNum);
    fissionMito.a = unifrnd(newMito(mitoNum).b,newMito(mitoNum).a/2);
    newMito(mitoNum).a = newMito(mitoNum).a-fissionMito.a;
    oldX = newMito(mitoNum).trueCoords(1)+(fissionMito.a*2)*cos(newMito(mitoNum).rotate);
    oldY = newMito(mitoNum).trueCoords(2)+(fissionMito.a*2)*sin(newMito(mitoNum).rotate);
    newX = fissionMito.trueCoords(1)-(travelDistanceAll(mitoNum)+newMito(mitoNum).a+fissionMito.a)*cos(newMito(mitoNum).rotate);
    newY = fissionMito.trueCoords(2)-(travelDistanceAll(mitoNum)+newMito(mitoNum).a+fissionMito.a)*sin(newMito(mitoNum).rotate);
    
    fissionMitoIm = zeros(frameSize,frameSize);

    %fill in the mitochondria and plot them.
    drawMito([newX,newY],normrnd(fissionMito.a*2,0.1),normrnd(fissionMito.b*2,0.1),fissionMito.rotate)

    %Make the plot into an image.
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    rect = [pos(1),pos(2),pos(3),pos(4)];
    f = getframe(gcf,rect);
    [X,~] = frame2im(f);
    
    %inverse colors for tracking.
    tempMitoIm = 1-X(:,:,1);
%     fissionMitoIm(:,:) = imgaussfilt(tempMitoIm*fissionMito.intensity,0.25/micronPerPixel);
    fissionMitoIm(:,:) = tempMitoIm*normrnd(fissionMito.intensity,0.1);
    hold off
    imshow(im)
    axis image off
    hold on
    
    %fill in the mitochondria and plot them.
    drawMito([oldX,oldY],normrnd(newMito(mitoNum).a*2,0.1),normrnd(newMito(mitoNum).b*2,0.1),newMito(mitoNum).rotate)

    %Make the plot into an image.
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    rect = [pos(1),pos(2),pos(3),pos(4)];
    f = getframe(gcf,rect);
    [X,~] = frame2im(f);
    
    %inverse colors for tracking.
    tempMitoIm = 1-X(:,:,1);
%     newMitoIm(:,:,mitoNum) = imgaussfilt(tempMitoIm*newMito(mitoNum).intensity,0.25/micronPerPixel);
    newMitoIm(:,:,mitoNum) = tempMitoIm*normrnd(newMito(mitoNum).intensity,0.1);
    

    %get mito properties (will be slightly different than generated and 
    %plotted stats)
    CCNew = bwconncomp(fissionMitoIm(:,:));
    newMitoTemp = regionprops(CCNew,fissionMitoIm(:,:),'Area','Centroid','Eccentricity','MaxFeretProperties','Extrema','MeanIntensity','WeightedCentroid','MajorAxisLength','MinorAxisLength','Perimeter','PixelIdxList');
    
    CCOld = bwconncomp(newMitoIm(:,:,mitoNum));
    oldMitoTemp = regionprops(CCOld,newMitoIm(:,:,mitoNum),'Area','Centroid','Eccentricity','MaxFeretProperties','Extrema','MeanIntensity','WeightedCentroid','MajorAxisLength','MinorAxisLength','Perimeter','PixelIdxList');
    
    %Correct the angles of the mitochondria and assign NN, frame, and
    %initiliaze other parameters.
    newMitoTemp.MaxFeretAngle = angle0to180(180-newMitoTemp.MaxFeretAngle);
    newMitoTemp.NN =0;
    newMitoTemp.NNExtrema = 0;
    newMitoTemp.frame = frameNum;
    newMitoTemp.NPA = 0;
    newMitoTemp.confident = 1;
    newMitoTemp.lost = 0;
    newMitoTemp.fission = mitoNum;
    newMitoTemp.fusion = 0;
    newMitoTemp.a = fissionMito.a;
    newMitoTemp.b = fissionMito.b;
    newMitoTemp.rotate = fissionMito.rotate;
    newMitoTemp.intensity = fissionMito.intensity;
    newMitoTemp.trueCoords = [newX,newY];
    newMitoTemp.overlap1 = fissionMito.overlap1;

   
    numMitoTotal = length(newMito);
    newMito(numMitoTotal+1) = newMitoTemp;
    newMitoIm(:,:,numMitoTotal+1) = fissionMitoIm;
    
    
    %Correct the angles of the mitochondria and assign NN, frame, and
    %initiliaze other parameters.
    oldMitoTemp.MaxFeretAngle = angle0to180(180-oldMitoTemp.MaxFeretAngle);
    oldMitoTemp.NN = newMito(mitoNum).NN;
    oldMitoTemp.NNExtrema = newMito(mitoNum).NNExtrema;
    oldMitoTemp.frame = frameNum;
    oldMitoTemp.NPA = newMito(mitoNum).NPA;
    oldMitoTemp.confident = newMito(mitoNum).confident;
    oldMitoTemp.lost = newMito(mitoNum).lost;
    oldMitoTemp.fission = newMito(mitoNum).fission;
    oldMitoTemp.fusion = newMito(mitoNum).fusion;
    oldMitoTemp.a = newMito(mitoNum).a;
    oldMitoTemp.b = newMito(mitoNum).b;
    oldMitoTemp.rotate = newMito(mitoNum).rotate;
    oldMitoTemp.intensity = newMito(mitoNum).intensity;
    oldMitoTemp.trueCoords = [oldX,oldY];
    oldMitoTemp.overlap1 = newMito(mitoNum).overlap1;

    
    newMito(mitoNum) = oldMitoTemp;
    
    hold off
end

end