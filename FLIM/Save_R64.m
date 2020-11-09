function im = Save_R64(frequency)

%Default args
if nargin<1 %default to 80MHz
    frequency = 80000000;
end

omega = 2*pi*frequency;

%Im order
% 1 - Intensity
% 2 - Harmonic 1 Phase (degrees)
% 3 - Harmonic 1 Modulation
% 4 - Harmonic 2 Phase (degrees)
% 5 - Harmonic 2 Modulation

[Im,Names,Path] = Open_R64();
for imageNum = 1:size(Im,2) 
    %cut off the edges because they are weird.
    im(imageNum).intensity = Im{1,imageNum}(2:end-2,2:end-2,1);
    if isempty(min(im(imageNum).intensity(im(imageNum).intensity>0)))
        continue
    end
    im(imageNum).intensity = im(imageNum).intensity/min(im(imageNum).intensity(im(imageNum).intensity>0));
    
    
    im(imageNum).phase = Im{1,imageNum}(2:end-2,2:end-2,2);
    im(imageNum).mod = Im{1,imageNum}(2:end-2,2:end-2,3);
    im(imageNum).phase2 = Im{1,imageNum}(2:end-2,2:end-2,2);
    im(imageNum).mod2 = Im{1,imageNum}(2:end-2,2:end-2,3);
    
    im(imageNum).tauPhase = tand(Im{1,imageNum}(2:end-2,2:end-2,2))/(omega);
    im(imageNum).tauMod = real(sqrt(((1./(Im{1,imageNum}(2:end-2,2:end-2,3))).^2)-1)./(omega));
    
    im(imageNum).g = Im{1,imageNum}(2:end-2,2:end-2,3).*cosd(Im{1,imageNum}(2:end-2,2:end-2,2));
    im(imageNum).s = Im{1,imageNum}(2:end-2,2:end-2,3).*sind(Im{1,imageNum}(2:end-2,2:end-2,2));
    
    im(imageNum).g2 = Im{1,imageNum}(2:end-2,2:end-2,5).*cosd(Im{1,imageNum}(2:end-2,2:end-2,4));
    im(imageNum).s2 = Im{1,imageNum}(2:end-2,2:end-2,5).*sind(Im{1,imageNum}(2:end-2,2:end-2,4));
    
    im(imageNum).name = Names{imageNum}; 
    im(imageNum).path = Path;
end
end