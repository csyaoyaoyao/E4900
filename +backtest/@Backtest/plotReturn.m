function plotReturn(obj)
% plot return curve
retnStratsTT = obj.BacktestResult;
f=figure;
f.Position(3)= 2*f.Position(3);
for i=1:size(retnStratsTT.Variables, 2)
    plot(retnStratsTT.Time, retnStratsTT{:,i}); hold on;
end
hold off;

datetick('x','mm/dd//yy','keepticks');
xlabel('Time');
ylabel('Return');
title('Portfolio Returns');
leg=legend(retnStratsTT.Properties.VariableNames, 'Location','NorthEastOutside');
set(leg,'Interpreter', 'none')

end

