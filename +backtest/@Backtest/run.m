function ResultTover = run(obj, assetTT, signalTT, varargin)
% Note
% cash is considered with return as RiskFreeRate.
% cost is considered to grow with CostTimeValueRate.

if any(any(ismissing(assetTT),1),2)
    error('Missing data in the dataset.');
end
%disp(assetTT)
[nTimeSteps, nAssets] = size(assetTT.Variables);

ip = inputParser;
ip.StructExpand = true;
ip.CaseSensitive = false;
ip.addParameter('Start', 1, ...
    @(x) validateattributes(x, {'double'}, {'scalar', 'integer', '>', 0, '<=', nTimeSteps}));
ip.addParameter('End', nTimeSteps, ...
    @(x) validateattributes(x, {'double'}, {'scalar', 'integer', '>', 0, '<=', nTimeSteps}));
ip.parse(varargin{:});
result = ip.Results;

% Debug if "Start" input and "End" input are valid
startStep = result.Start;
endStep= result.End;
if startStep > endStep
    error('Starting step should be less than ending step.');
end

strategies = obj.Strategies;

backTestRange = startStep:endStep;
nStrategies = numel(strategies);
nBackTestSteps = numel(backTestRange);

% extract strategy names
stratsName = arrayfun(@(x)x.Name, strategies, 'UniformOutput',false);

% Initiate all variables to store intermediary results
% 

obj.InitialPosition = signalTT{1,:}.';

retnSeries = zeros(nBackTestSteps, nStrategies);
turnover = zeros(nBackTestSteps, nStrategies);
buycost =  zeros(nBackTestSteps, nStrategies);
sellcost =  zeros(nBackTestSteps, nStrategies);
managementCost =  zeros(nBackTestSteps, nStrategies);
performanceCost =  zeros(nBackTestSteps, nStrategies);
pwgts = repmat(obj.InitialPosition, 1, nStrategies);  % weights from last rebalancing step
highwatermark =  ones(1, nStrategies)*obj.InitialMarketValue;
endAcctVal = ones(1, nStrategies)*obj.InitialMarketValue;

% Iterate over all strategies
% Initiate the intermediary results and output variables
for i = 1:numel(strategies)
    stratName = strategies(i).Name;
    pwgtStructSOD.(stratName) = zeros(nBackTestSteps, nAssets+1);
    positionStructSOD.(stratName) = zeros(nBackTestSteps, nAssets+1);
    pwgtStructEOD.(stratName) = zeros(nBackTestSteps, nAssets+1);
    positionStructEOD.(stratName) = zeros(nBackTestSteps, nAssets+1);
end
% previous market value for cash + asset positions. [cash; asset_pwgts] (1+N)-by-1 vector
prevMktPosi = obj.InitialMarketValue*repmat([1-sum(obj.InitialPosition); obj.InitialPosition], ...
    1, nStrategies);
cashIdx = 1;
assetIdxs = 2:(numel(obj.InitialPosition)+1);
% main logic for backtesting
for iStep = 2:nBackTestSteps
    idx = backTestRange(iStep);
    assetCurrData = assetTT{idx, :};
    for i = 1:numel(strategies)
        currStrat = strategies(i);
        sodMktPosi = prevMktPosi(:, i);
        sodPwgt = sodMktPosi(assetIdxs)./sum(sodMktPosi);   % sod actual pwgt given previous market change
        sodMktVal = sum(sodMktPosi);       
        if currStrat.RebalanceFreq==1 || mod(iStep, currStrat.RebalanceFreq) == 1
            if idx-currStrat.LookbackStep>0
                assetLookbackData = assetTT{idx-currStrat.LookbackStep:idx-1, :};
                signalLookbackData = signalTT{idx-currStrat.LookbackStep:idx-1, :};
                sodPwgt_ = sodPwgt';
                sodPwgtUpdated = currStrat.AllocationFcnHandle(assetLookbackData, signalLookbackData, sodPwgt_);
                sodPwgtUpdated = sodPwgtUpdated/sum(sodPwgtUpdated);
                % pwgt = run(obj, retn, signal, currPwgt)
                turnover(iStep, i) = sum(abs(sodPwgtUpdated - sodPwgt));  % both buy & sell with wgts to compute Turnover
                [buycost(iStep, i) , sellcost(iStep, i)] = currStrat.computeTransactionCost(sodMktVal*(sodPwgtUpdated - sodPwgt));
                pwgts(:, i) = sodPwgtUpdated;
            end
        else
            pwgts(:, i) = sodPwgt;
        end
        
        retn = [obj.RiskFreeRate assetCurrData];
        % SOD : Start of day
        % update SOD mkt value for each position given new pwgts
        cash_assets_pwgt = [1-sum(pwgts(:, i)); pwgts(:, i)];
        sodMktPosi = sodMktVal*cash_assets_pwgt';
        pwgtStructSOD.(strategies(i).Name)(iStep, :) = cash_assets_pwgt';
        positionStructSOD.(strategies(i).Name)(iStep, :) = sodMktPosi;
        
        % EOD : End of Day
        % udpate EOD mkt value for each position with current day's market change
        eodMktPosi = sodMktPosi.*(1+retn);
        pwgtStructEOD.(strategies(i).Name)(iStep, :) = eodMktPosi/sum(eodMktPosi);
        positionStructEOD.(strategies(i).Name)(iStep, :) = eodMktPosi;
        
        % keep all costs apart, and account for them when computing the returns
        % managementCost, performanceCost, buycost, sellcost store the total costs from step 1 to iStep
        if currStrat.ManagementFeeFreq==1 || mod(iStep, currStrat.ManagementFeeFreq) == 1           
            eodMktVal = sum(eodMktPosi);
            endAcctVal(i) = eodMktVal;
            managementCost(iStep, i) = currStrat.computeManagementCost(endAcctVal(i));
            performanceCost(iStep, i) = currStrat.computePerformanceCost(highwatermark(i), eodMktVal);
            highwatermark(i) = max(highwatermark(i), eodMktVal);
            endAcctVal(i) = eodMktVal;
        end
        
        buycost(iStep, i) = buycost(iStep, i)*(1+obj.CostTimeValueRate);
        sellcost(iStep, i) = sellcost(iStep, i)*(1+obj.CostTimeValueRate);
        if iStep>1
            buycost(iStep, i) = buycost(iStep, i) + buycost(iStep-1, i)*(1+obj.CostTimeValueRate);
            sellcost(iStep, i) = sellcost(iStep, i) + sellcost(iStep-1, i)*(1+obj.CostTimeValueRate);
            
            managementCost(iStep, i) = managementCost(iStep, i)+ managementCost(iStep-1, i)*(1+obj.CostTimeValueRate);  
            performanceCost(iStep, i) = performanceCost(iStep, i) + performanceCost(iStep-1, i)*(1+obj.CostTimeValueRate);  
        end
        
        % decrease SOD and EOD mktval by cost
        if iStep>1
            sodMktVal = sodMktVal - buycost(iStep-1, i)- sellcost(iStep-1, i)-...
                managementCost(iStep-1, i) - performanceCost(iStep-1, i);
        end
        eodMktVal = sum(eodMktPosi) - buycost(iStep, i) - sellcost(iStep, i)-...
            managementCost(iStep, i) - performanceCost(iStep, i);
        
        retnSeries(iStep, i) = eodMktVal/sodMktVal-1;
        prevMktPosi(:, i) = eodMktPosi;
    end
end

for i = 1:numel(strategies)
    stratName = strategies(i).Name;
    pwgtSeries = array2table([pwgtStructSOD.(stratName)],...
        'VariableNames',["Cash", assetTT.Properties.VariableNames]);
    Time = assetTT.Time([backTestRange]);
    pwgtSeries = addvars(pwgtSeries, Time, 'Before', "Cash");
    obj.PwgtStructSOD.(stratName) = table2timetable(pwgtSeries);
    
    posiSeries = array2table([positionStructSOD.(stratName)],...
        'VariableNames',["Cash", assetTT.Properties.VariableNames]);
    Time = assetTT.Time([backTestRange]);
    posiSeries = addvars(posiSeries, Time, 'Before', "Cash");
    obj.PositionStructSOD.(stratName) = table2timetable(posiSeries);
    
    pwgtSeries = array2table([pwgtStructEOD.(stratName)],...
        'VariableNames',["Cash", assetTT.Properties.VariableNames]);
    Time = assetTT.Time([backTestRange]);
    pwgtSeries = addvars(pwgtSeries, Time, 'Before', "Cash");
    obj.PwgtStructEOD.(stratName) = table2timetable(pwgtSeries);
    
    posiSeries = array2table([positionStructEOD.(stratName)],...
        'VariableNames',["Cash", assetTT.Properties.VariableNames]);
    Time = assetTT.Time([backTestRange]);
    posiSeries = addvars(posiSeries, Time, 'Before', "Cash");
    obj.PositionStructEOD.(stratName) = table2timetable(posiSeries);
end

retnSeries = array2table([retnSeries],...
    'VariableNames',stratsName);
Time = assetTT.Time([backTestRange]);
retnSeries = addvars(retnSeries, Time,'Before',stratsName{1});
ResultTT = table2timetable(retnSeries);
obj.BacktestResult = ResultTT;

turnover = array2table([turnover],...
    'VariableNames',stratsName);
Time = assetTT.Time([backTestRange]);
turnover = addvars(turnover, Time,'Before',stratsName{1});
ResultTover = table2timetable(turnover);
obj.Turnover = ResultTover ; %table2timetable(turnover);

buycost = array2table([buycost],...
    'VariableNames',stratsName);
Time = assetTT.Time([backTestRange]);
buycost = addvars(buycost, Time,'Before',stratsName{1});
obj.BuyCost = table2timetable(buycost);

sellcost = array2table([sellcost],...
    'VariableNames',stratsName);
Time = assetTT.Time([backTestRange]);
sellcost = addvars(sellcost, Time,'Before',stratsName{1});
obj.SellCost = table2timetable(sellcost);

managementCost = array2table([managementCost],...
    'VariableNames',stratsName);
Time = assetTT.Time([backTestRange]);
managementCost = addvars(managementCost, Time,'Before',stratsName{1});
obj.ManagementCost = table2timetable(managementCost);

performanceCost = array2table([performanceCost],...
    'VariableNames',stratsName);
Time = assetTT.Time([backTestRange]);
performanceCost = addvars(performanceCost, Time,'Before',stratsName{1});
obj.PerformanceCost = table2timetable(performanceCost);
end