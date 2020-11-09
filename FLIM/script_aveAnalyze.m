allNADH = Save_R64();
allTMRM = Save_R64();
%%

noiseThreshold = 0.25;

highThresh = 30;
lowThresh = 15;

int = allNADH.intensity;
int(allNADH.intensity>highThresh) = nan;
int(allNADH.intensity<lowThresh) = nan;
cyto = int;
cytog = allNADH.g;
cytog(isnan(cyto)) = nan;
cytos = allNADH.s;
cytos(isnan(cyto)) = nan;

int = allNADH.intensity;
int(allNADH.intensity<highThresh) = nan;
mito = int;
mitog = allNADH.g;
mitog(allNADH.intensity<highThresh) = nan;
mitos = allNADH.s;
mitos(allNADH.intensity<highThresh) = nan;

int = allNADH.intensity;
int(allNADH.intensity<lowThresh) = nan;
all = int;
allg = allNADH.g;
allg(allNADH.intensity<lowThresh) = nan;
alls = allNADH.s;
alls(allNADH.intensity<lowThresh) = nan;

cytofb = [];

[frac1,frac2,frac3] = fractionCalc(cytos(cytog>0),cytog(cytog>0),0.4e-09,3.2e-09,80000000);
% frac1 = frac1(abs(frac3)<noiseThreshold);
% frac2 = frac2(abs(frac3)<noiseThreshold);
% frac3 = frac3(abs(frac3)<noiseThreshold);
cytofb(:,1) = frac1./(frac1+frac2);
cytofb(:,2) = frac2./(frac1+frac2);
cytofb(:,3) = frac3;

mitofb = [];

[frac1,frac2,frac3] = fractionCalc(mitos(mitog>0),mitog(mitog>0),0.4e-09,3.2e-09,80000000);
% frac1 = frac1(abs(frac3)<noiseThreshold);
% frac2 = frac2(abs(frac3)<noiseThreshold);
% frac3 = frac3(abs(frac3)<noiseThreshold);
mitofb(:,1) = frac1./(frac1+frac2);
mitofb(:,2) = frac2./(frac1+frac2);
mitofb(:,3) = frac3;

allfb = [];

[frac1,frac2,frac3] = fractionCalc(alls(allg>0),allg(allg>0),0.4e-09,3.2e-09,80000000);
% frac1 = frac1(abs(frac3)<noiseThreshold);
% frac2 = frac2(abs(frac3)<noiseThreshold);
% frac3 = frac3(abs(frac3)<noiseThreshold);
allfb(:,1) = frac1./(frac1+frac2);
allfb(:,2) = frac2./(frac1+frac2);
allfb(:,3) = frac3;