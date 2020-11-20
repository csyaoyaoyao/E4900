function plotTransactionCost(obj)
% plot return curve
buycost = obj.BuyCost;
sellCost = obj.SellCost;
f=figure;
f.Position(3)= 2*f.Position(3);
for i=1:size(buycost.Variables, 2)
    bar(buycost.Time, buycost{:,i}); hold on;
    bar(sellCost.Time, -sellCost{:,i}); hold on;
end
hold off;

datetick('x','mm/dd//yy','keepticks');
xlabel('Time');
ylabel('Return');
title('Portfolio Returns');
leg=legend(buycost.Properties.VariableNames, 'Location','NorthEastOutside');
set(leg,'Interpreter', 'none')

end

