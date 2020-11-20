function pwgt= allocByEqualWeight(data, signal, currPwgt)
% equal weighted allocation
NAssets = size(data, 2);
if signal==0
   pwgt = zeros(NAssets, 1);
   return
end
pwgt = ones(NAssets, 1);
pwgt = pwgt/sum(pwgt);
end