fileID = fopen('/Users/aelefebv/Desktop/3D/fusion.txt','r');
fusion= fscanf(fileID,'%f');
%%
close
histogram(sol(A==0),200,'Normalization','cdf')
hold on
histogram(sol(A==1),200,'Normalization','cdf')
histogram(sol(A==2),200,'Normalization','cdf')
%%
close
histogram((majAx(A==0)./minAx(A==0)),50,'Normalization','cdf')
hold on
histogram((majAx(A==1)./minAx(A==1)),50,'Normalization','cdf')
histogram((majAx(A==2)./minAx(A==2)),50,'Normalization','cdf')
%%
close
histogram(int(A==0),200,'Normalization','pdf')
hold on
histogram(int(A==1),200,'Normalization','pdf')
% histogram(int(A==2),200,'Normalization','pdf')
%%
var = log10(vol);
close
histogram(var(A==0),100,'Normalization','probability')
hold on
histogram(var(A==1),100,'Normalization','probability')
histogram(var(A==2),100,'Normalization','probability')
%%
var = meanSpeed;
close
histogram(var(A==0),400,'Normalization','cdf')
hold on
histogram(var(A==1),400,'Normalization','cdf')
histogram(var(A==2),400,'Normalization','cdf')
%%
var = meanSpeed;
close
histogram(log10(var(A==0)),400,'Normalization','pdf')
hold on
histogram(log10(var(A==1)),400,'Normalization','pdf')
histogram(log10(var(A==2)),400,'Normalization','pdf')

getCoV(var(A==0))
getCoV(var(A==1))
getCoV(var(A==2))

%%

numFisLow = nnz(fission(A==0));
numFisMed = nnz(fission(A==1));
numFisHigh = nnz(fission(A==2));

numTotLow = numel(fission(A==0));
numTotMed = numel(fission(A==1));
numTotHigh = numel(fission(A==2));

fracFisLow = numFisLow/numTotLow;
fracFisMed = numFisMed/numTotMed;
fracFisHigh = numFisHigh/numTotHigh;
%%

numFusLow = nnz(fusion(A==0));
numFusMed = nnz(fusion(A==1));
numFusHigh = nnz(fusion(A==2));

numTotLow = numel(fusion(A==0));
numTotMed = numel(fusion(A==1));
numTotHigh = numel(fusion(A==2));

fracFusLow = numFusLow/numTotLow;
fracFusMed = numFusMed/numTotMed;
fracFusHigh = numFusHigh/numTotHigh;


%%
allOut = [majAx./minAx,sol,peri,meanSpeed,dir,int,A];
allOutHeaders = {'AxRatio', 'Sol', 'Peri', 'MeanSpeed', 'Dir', 'Int','Invasivity'};
writeOut = [allOutHeaders;num2cell(allOut)];
writecell(writeOut,'/Users/aelefebv/Desktop/2D/params2D.txt');
