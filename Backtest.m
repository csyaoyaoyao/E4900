classdef Backtest < handle
    %BACKTEST takes the input strategies and apply them on a given time
    % period to allow the visualization and comparison of the overall 
    % performances of the strategies.
    
    % Notes on enhancements:
    % -- Input validation, warning and error messages.
    % -- Rolling window approach is used to apply allocation strategies.
    % Potentially we can also support walk forward approach with expanding
    % windows.
    % -- The costs can either not grow or grow with a certain rate. This is
    % controled by CostTimeValueRate.    
    % -- May need to track the value of highwater mark and offer its visualization.
    
    properties
        Strategies
        RiskFreeRate
        CostTimeValueRate
        InitialMarketValue
        InitialPosition
        NumAssets
    end
       
    properties(SetAccess=private)
        BacktestResult  % a timetable with returns from different strategies
        PwgtStructSOD
        PositionStructSOD
        PwgtStructEOD
        PositionStructEOD        
        Turnover        % a timetable with turnovers from different strategies
        BuyCost
        SellCost
        ManagementCost  % proportional to the portfolio value when evaluating the cost
        PerformanceCost % proportional to the increase compared to highwater mark
    end
    
    methods
        function obj = Backtest(strategies, numAssets, varargin)
            obj.Strategies = strategies;
            obj.NumAssets = numAssets;
            
            ip = inputParser;
            ip.StructExpand = true;
            ip.CaseSensitive = false;
            ip.addParameter('RiskFreeRate', 0, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', 'real'}));  
            ip.addParameter('CostTimeValueRate', 0, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', 'real'}));              
            ip.addParameter('InitialMarketValue', 1e4, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', 'real'}));        
            % Initial position: defaulted to be zero position
            % cashPosition = 1-sum(initialPosition)
            ip.addParameter('InitialPosition', zeros(obj.NumAssets, 1)/obj.NumAssets, ...
                @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', obj.NumAssets})); 
            
            ip.parse(varargin{:});
            result = ip.Results;
            obj.InitialPosition = result.InitialPosition;       % the sum can be less than 1      
            obj.RiskFreeRate = result.RiskFreeRate;  
            obj.CostTimeValueRate = result.CostTimeValueRate;
            obj.InitialMarketValue = result.InitialMarketValue;
          
        end  
        
        plotReturn(obj);
        plotEquityCurve(obj, varargin);
        plotTransactionCost(obj);
        ResultTT = run(obj, assetTT, signalTT, varargin);
        PerfTT = summary(obj);
        printSummary(obj);
    end
end

