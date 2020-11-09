function newFusionMito = simulateFusion(currentMito,fusionMito)
newFusionMito = fusionMito;
newFusionMito.a = currentMito.a + fusionMito.a;
newFusionMito.WeightedCentroid(1) = (currentMito.WeightedCentroid(1)*currentMito.a + fusionMito.WeightedCentroid(1)*fusionMito.a)/(currentMito.a+fusionMito.a);
newFusionMito.WeightedCentroid(2) = (currentMito.WeightedCentroid(2)*currentMito.a + fusionMito.WeightedCentroid(2)*fusionMito.a)/(currentMito.a+fusionMito.a);
end