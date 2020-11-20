function pwgt= allocByMeanVariancePortOptim(dataTT, signal, currPwgt)
numAssets = size(dataTT, 2);
if signal==0
   pwgt = zeros(numAssets, 1);
   return
end
p = Portfolio('NumAssets', numAssets, 'Budget', 1, 'LowerBound', 0, ...
    'UpperBound', 1);
p = estimateAssetMoments(p, dataTT);
targetRetn = mean(p.estimatePortReturn(p.estimateFrontierLimits('both')));
pwgt = p.estimateFrontierByReturn(targetRetn);

end