function newFusionMito = simulateFusion(currentMito,fusionMito)
newFusionMito = fusionMito;
newFusionMito.a = currentMito.a + fusionMito.a;
newFusionMito.WeightedCentroid(1) = mean(currentMito.WeightedCentroid(1) + fusionMito.WeightedCentroid(1));
newFusionMito.WeightedCentroid(2) = mean(currentMito.WeightedCentroid(2) + fusionMito.WeightedCentroid(2));
end