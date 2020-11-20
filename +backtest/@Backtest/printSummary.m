function printSummary(obj)
result = obj.BacktestResult;
fprintf(['~~~~~~~~~~~~ Backtest Summary ~~~~~~~~~~~~ \n']);

fprintf('\n');
fprintf('Start                  : %s\n', result.Time(1));
fprintf('End                    : %s\n', result.Time(end));
fprintf('Number of Strategies   : %d\n\n', numel(result.Properties.VariableNames));

backTestRetns = result.Variables;
MeanRetn = mean(backTestRetns);
StdRetn = sqrt(var(backTestRetns));
SharpeRatio = MeanRetn./StdRetn;

cmpRetn = ret2tick(backTestRetns, 'method', 'simple');
MaxDD = maxdrawdown(cmpRetn);
TotalRetn = cmpRetn(end, :)-1;
avgTO = mean(obj.Turnover.Variables);
avgBuyCost = mean(obj.BuyCost.Variables);  % averaged over trading period
avgSellCost = mean(obj.SellCost.Variables);
managementCost = obj.ManagementCost{end, :};

strats = pad(result.Properties.VariableNames);

fprintf('-- Total Return (%%)--       \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, TotalRetn(i)*100);
end
fprintf('\n');

fprintf('-- Average Turnover (%%)--      \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, avgTO(i)*100);
end
fprintf('\n');

fprintf('-- Average BuyCost (Dollar Amt)--      \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, avgBuyCost(i));
end
fprintf('\n');

fprintf('-- Average SellCost (Dollar Amt)--      \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, avgSellCost(i));
end
fprintf('\n');

fprintf('-- Sharpe Ratio --       \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, SharpeRatio(i));
end
fprintf('\n');
fprintf('-- Average Return (%%)--    \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, MeanRetn(i)*100);
end
fprintf('\n');

fprintf('-- Volatility (%%)--         \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, StdRetn(i)*100);
end
fprintf('\n');

fprintf('-- Max Drawdown (%%)--      \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, MaxDD(i)*100);
end
fprintf('\n');

fprintf('-- Management Fee (Dollar Amt)--      \n');
for i = 1:numel(strats)
    fprintf('%s : %8.4f\n', strats{i}, managementCost(i));
end
fprintf('\n');