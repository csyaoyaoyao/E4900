function out = riskCostFunc(wts, covMat)
    nAssets = numel(wts);
    portVar = wts*covMat*wts';
    margRisk = wts' .* (covMat * wts') ./ sqrt(portVar);
    dif = margRisk - sqrt(portVar)/nAssets;  % equal marginal contribution to risk
    out = 0.5*sum(dif.^2);
end