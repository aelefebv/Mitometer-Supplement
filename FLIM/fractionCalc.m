function [fraction1,fraction2,fraction3] = fractionCalc(s,g,freeTau,boundTau,frequency)

omega = 2*pi*frequency;

g1 = 1/(1+(omega*freeTau)^2);
s1 = (omega*freeTau)/(1+(omega*freeTau)^2);

g2 = 1/(1+(omega*boundTau)^2);
s2 = (omega*boundTau)/(1+(omega*boundTau)^2);

A =         [g1  g2  0;
             s1  s2  0;
             1      1       1];
         
x = [];
for pixel = 1:numel(g)
    y = [g(pixel);s(pixel);1];
    x(:,pixel) = A\y;
end

pixelNum = 1;
for j = 1:numel(g)
    fraction1(j) = x(1,pixelNum);
    fraction2(j) = x(2,pixelNum);
    fraction3(j) = x(3,pixelNum);
    pixelNum = pixelNum+1;
end
if ~numel(g)
    fraction1 = nan;
    fraction2 = nan;
    fraction3 = nan;
end

end