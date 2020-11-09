function CoV = getCoV(var)

meanVar = nanmean(var);
stdVar = nanstd(var);
CoV = stdVar/meanVar;
end