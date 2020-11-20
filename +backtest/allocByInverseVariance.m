function pwgt= allocByInverseVariance(data, signal, currPwgt)
% inverse-variance allocation
if signal==0
   pwgt = zeros(size(data, 2), 1);
   return
end
assetCov = cov(data);
pwgt = 1./diag(assetCov);
pwgt = pwgt/sum(pwgt);
end