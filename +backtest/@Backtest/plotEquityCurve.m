function plotEquityCurve(obj, varargin)
% plot equity curve
retnStratsTT = obj.BacktestResult;

ip = inputParser;
ip.StructExpand = true;
ip.CaseSensitive = false;
ip.addParameter('IncludeCost', true, ...
    @(x) validateattributes(x, {'logical'}, {'scalar'}));
ip.parse(varargin{:});
result = ip.Results;

includeCost = result.IncludeCost;

% retn->tick with simple method
equityData = ret2tick(retnStratsTT, 'StartPrice', obj.InitialMarketValue, 'method', 'simple');
f=figure;
f.Position(3)= 2*f.Position(3);
for i=1:size(retnStratsTT.Variables, 2)
    if includeCost
        plot(equityData.Time, equityData{:,i}); hold on;
    else
        mvTT = obj.PositionStructSOD.(obj.Strategies(i).Name);
        mktval = sum(mvTT.Variables, 2);
        plot(mvTT.Time, mktval); hold on;        
    end
end
    
hold off;

datetick('x','mm/dd//yy','keepticks');
xlabel('Time');
ylabel('Value');
title('Portfolio Equity Curve');
leg=legend(arrayfun(@(x) x.Name, obj.Strategies, 'UniformOutput', false),...
    'Location','NorthEastOutside');
set(leg,'Interpreter', 'none')

end

