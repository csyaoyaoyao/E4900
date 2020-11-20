function [pwgt, vari] = allocByEqualRiskContribution(dataTT, signal, currPwgt)
% perform equal markginal risk contribution analysis
covRet = cov(dataTT);
nAsset = size(covRet, 2);
if signal==0
   pwgt = zeros(nAsset, 1);
   vari = 0;
   return
end
% Set up the function handle
ercObj = @(x) riskCostFunc(x , covRet);
options = optimoptions('fmincon','Algorithm','interior-point', 'Display','off', 'OptimalityTolerance', 1e-8);
x0 = ones(1 , nAsset)./ nAsset;
Aeq = ones(1 , nAsset);
Beq = 1;
LB = zeros(nAsset, 1);
UB = ones(nAsset, 1);
pwgt = fmincon(ercObj, x0, [], [], Aeq, Beq, LB, UB, [], options);
pwgt = pwgt';
vari = pwgt'*covRet*pwgt;
    function out = riskCostFunc(wts, covMat)
        nAssets = numel(wts);
        portVar = wts*covMat*wts';
        margRisk = wts' .* (covMat * wts') ./ sqrt(portVar);
        dif = margRisk - sqrt(portVar)/nAssets;  % equal marginal contribution to risk
        out = 0.5*sum(dif.^2);
    end
end