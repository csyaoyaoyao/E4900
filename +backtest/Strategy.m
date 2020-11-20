classdef Strategy
    % Make decisions on the asset allocations in the investment universe.
    % Maintain the states related to the strategy.
    % Store properties specific to the strategy
    
    properties
        Name
        LookbackStep
        RebalanceFreq
        AllocationFcnHandle        % logic to decide allocations among assets
        
        % optional inputs
        TransactionCostFcnHandle   % logic to decide the buy/sell cost given the change of the assets' market value
        MarketImpactFcnHandle
        ManagementFeeRate
        PerformanceFeeRate
        ManagementFeeFreq
    end
    
    % maintain the states
    properties
       CurrState 
    end
    
    methods(Static)
        function [buy, sell] = computeTransactionCostDefault(deltaMktValue)
            numAssets = numel(deltaMktValue);
            buy = sum(max(deltaMktValue, zeros(numAssets,1)))*0.01;
            sell = sum(max(0-deltaMktValue, zeros(numAssets,1)))*0.01;
        end
        
        function impactCost = computeMarketImpactDefault(deltaMktValue)
            impactCost = 0*deltaMktValue;
        end
    end
    
    methods
        function obj = Strategy(name, lookbackStep, rebalanceFreq, ...
                fcnHandle, varargin)
            obj.Name = name;
            obj.LookbackStep = lookbackStep;
            obj.RebalanceFreq = rebalanceFreq; 
            obj.AllocationFcnHandle = fcnHandle;
            
            % flexible in terms of allowing a custom cost analyzer
            % function, and management fee, etc. 
            defaultCostFcnHandle = @(x) backtest.Strategy.computeTransactionCostDefault(x);
            defaultMarketImpactFcnHandle = @(x) backtest.Strategy.computeMarketImpactDefault(x);
            ip = inputParser;
            ip.StructExpand = true;
            ip.CaseSensitive = false;
            ip.addParameter('TransactionCostFcnHandle', defaultCostFcnHandle);  
            ip.addParameter('MarketImpactFcnHandle', defaultMarketImpactFcnHandle); 
            ip.addParameter('ManagementFeeRate', 0, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', 'real'}));    
            ip.addParameter('PerformanceFeeRate', 0, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', 'real'}));  
            ip.addParameter('ManagementFeeFreq', 1, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', 'finite', 'real'})); 
            ip.parse(varargin{:});
            result = ip.Results;
            obj.TransactionCostFcnHandle = result.TransactionCostFcnHandle; 
            obj.ManagementFeeRate = result.ManagementFeeRate;
            obj.PerformanceFeeRate = result.PerformanceFeeRate; 
            obj.ManagementFeeFreq = result.ManagementFeeFreq;
        end
        
        function pwgt = run(obj, retn, signal, currPwgt)
            % allocate among assets
            % retn: a matrix with asset returns
            % signal: a matrix with investment signals
            % currPwgt: a vector with asset weights
            pwgt = obj.AllocationFcnHandle(retn, signal, currPwgt);
        end
        
        function [buy, sell] = computeTransactionCost(obj, deltaMktValue)
            [buy, sell] = obj.TransactionCostFcnHandle(deltaMktValue);
        end
        
        function performanceCost = computePerformanceCost(obj, highwatermark, currPortVal)
            % highwatermark: a scalar that records the highest portfolio
            % value during the backtesting.
            %
            % currPortVal: a vector of positions (in market value), instead
            % of directly passing in the total mkt value of the portfolio.
            % Because a user may define their performance cost that is
            % asset specific. 
            
            performanceCost = 0;
            currPortTotalVal = sum(currPortVal);
            if  currPortTotalVal> highwatermark
                performanceCost = (currPortTotalVal - highwatermark)*obj.PerformanceFeeRate;
            end
        end
        
        function managementCost = computeManagementCost(obj, endPortfVal)
            % endPortfVal: a vector of positions at the end of current trading 
            % (in market value), instead of directly passing in the total mkt value 
            % of the portfolio. Because a user may define their performance cost
            % that is asset specific. 
            managementCost = sum(endPortfVal)*obj.ManagementFeeRate*obj.ManagementFeeFreq;
        end

    end
end

