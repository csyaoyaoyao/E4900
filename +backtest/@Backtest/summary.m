function performanceTable = summary(obj)
% compute performance metric
backTestRetns = obj.BacktestResult.Variables;
MeanRetn = mean(backTestRetns);
StdRetn = sqrt(var(backTestRetns));
SharpeRatio = (MeanRetn-obj.RiskFreeRate)./StdRetn;

cmpRetn = ret2tick(backTestRetns, 'method', 'simple');
MaxDD = maxdrawdown(cmpRetn);
TotalRetn = cmpRetn(end, :)-1;
avgTO = mean(obj.Turnover.Variables);
avgBuyCost = mean(obj.BuyCost.Variables);
avgSellCost = mean(obj.SellCost.Variables);
maxTO = max(obj.Turnover.Variables);
managementCost = obj.ManagementCost{end, :};
performanceCost = obj.PerformanceCost{end, :};
perfData = [TotalRetn', avgTO', maxTO', SharpeRatio', MeanRetn', StdRetn', ...
    MaxDD', avgBuyCost', avgSellCost', managementCost', performanceCost']';

metricNames = {'TotalRetn', 'AverageTurnover', ...
    'MaxTurnover', 'SharpeRatio', 'AverageRetn', 'Volatility', 'MaxDrawdown', ...
    'AverageBuyCost', 'AverageSellCost', 'ManagementFee', 'PerformanceCost'};
colStratsNames = (obj.BacktestResult.Properties.VariableNames)';


performanceTable = array2table(perfData, 'VariableNames', colStratsNames, 'RowNames', metricNames);

end