classdef backtestStrategy <  matlab.mixin.Heterogeneous
    %BACKTESTSTRATEGY Strategy for portfolio backtesting.
    %
    % Syntax:
    %
    %   strat = backtestStrategy(name, rebalanceFcn)
    %   strat = backtestStrategy(name, rebalanceFcn, param1, value1,...)
    %
    % Description:
    %
    %   BACKTESTSTRATEGY takes as input a name and a function handle which
    %   specifies how the strategy will rebalance a portfolio of assets.
    %   BACKTESTSTRATEGY objects are used with the BACKTESTENGINE class to
    %   backtest portfolio trading strategies against historical market
    %   data.
    %
    % Input arguments:
    %
    %   name - Strategy name.  A string identifying the strategy.  Cannot
    %     contain special characters except for underscores.
    %
    %   rebalanceFcn - A function handle that computes new portfolio
    %     weights during the backtest.  The rebalanceFcn implements the
    %     core logic of the trading strategy and must have one of the
    %     following signatures:
    %
    %     new_weights = rebalanceFcn(weights,assetPrices)
    %     new_weights = rebalanceFcn(weights,assetPrices,signalData)
    %
    %     The rebalance function is called by the backtestEngine each time
    %     the strategy must be rebalanced.  The backtestEngine calls the
    %     rebalance function with the following arguments:
    %
    %       * weights: the current portfolio weights before rebalancing,
    %         specified as decimal percentages.
    %       * assetPrices: a timetable containing a rolling window of
    %         adjusted asset prices.
    %       * signalData: a timetable containing a rolling window of signal
    %         data.  If signal data is provided to the backtestEngine then
    %         the engine object will pass it to the strategy rebalance
    %         function (3 input argument syntax).  If signal data is not
    %         provided to the backtestEngine then the rebalance function
    %         will be called with the 2 input argument syntax.
    %
    %     The rebalance function must return a single output argument:
    %
    %       * new_weights: a vector of asset weights specified as decimal
    %         percentages.  If the new weights sum to 1 then the portfolio
    %         is fully invested.  If the weights sum to less than 1 then
    %         the portfolio will have the remainder in cash, earning the
    %         RiskFreeRate (specified in the backtestEngine).  If the
    %         weights sum to more than 1 then there will be a negative cash
    %         position (margin) and the cash borrowed will accrue interest
    %         at the CashBorrowRate (specified in the backtestEngine).
    %
    % Optional Input Parameter Name/Value Pairs:
    %
    %   'RebalanceFrequency' : A scalar that specifies how often the
    %     portfolio is rebalanced during the backtest.  The rebalance
    %     frequency is specified in terms of the number of time steps
    %     between rebalancing.  For example if the backtestEngine is
    %     provided with daily price data, then the RebalanceFrequency would
    %     specify the number of days between rebalancing.  The default is
    %     1, meaning the strategy would rebalance with each time step.
    %
    %   'TransactionCosts' : The transaction costs for trades.  Can be
    %     specified in the following ways:
    %
    %       * rate: A scalar decimal percentage charge to both purchases
    %           and sales of assets.  For example if TransactionCosts were
    %           set to 0.001, then each transaction (buys and sells) would
    %           pay 0.1% in transaction fees.
    %       * [buyRate, sellRate]: A 1-by-2 vector of decimal percentage
    %           rates that specifies separate rates for buying vs. selling
    %           of assets.
    %       * computeTransactionCostsFcn: A function handle to compute
    %           customized transaction costs.  If specified as a function
    %           handle, the backtestEngine will call the TransactionCosts
    %           function to compute the fees for each rebalance.  The
    %           provided function must have the following signature:
    %
    %             [buyCosts,sellCosts] = computeCostsFcn(deltaPositions)
    %
    %           It must take a single input argument, deltaPositions, which
    %           is a vector of changes in asset positions for all assets
    %           (in currency units) as a result of a rebalance.  Positive
    %           elements in the deltaPositions vector indicate purchases
    %           while negative entries represent sales.  The function
    %           handle must return two output arguments, buyCosts and
    %           sellCosts, which contain the total costs (in currency) for
    %           the entire rebalance for each type of transaction.
    %
    %     The default TransactionCosts is 0, mean transaction costs are not
    %     computed.
    %
    %   'LookbackWindow' : A 1-by-2 vector that specifies the minimum and
    %       maximum size of the rolling windows of data (asset prices and
    %       signal data) that are provided to the rebalance function.
    %       These limits are specified in terms of the number of time
    %       steps.  For example if the backtestEngine is provided with
    %       daily price data, then the lookback window would specify the
    %       size bounds of the rolling window in days.  The default is [0
    %       Inf], meaning all available past data is given to the rebalance
    %       function.
    %
    %       If a non-zero minimum is specified, then the rebalance function
    %       will not be called until enough time steps have processed to
    %       meet the minimum size.
    %
    %       Alternatively, the LookbackWindow can be set to a single scalar
    %       value indicating that the rolling window should be exactly that
    %       size.  The minimum and maximum size will both be set to the
    %       provided value.
    %
    %   'InitialWeights' : A vector of initial portfolio weights.  The
    %       InitialWeights vector sets the portfolio weights before the
    %       backtestEngine begins the backtest.  The size of the initial
    %       weights vector must match the number of assets used in the
    %       backtest.
    %
    %       Alternatively, you can set the InitialWeights to empty ([]),
    %       indicating the strategy will begin uninvested and in a 100%
    %       cash position.  The default is empty ([]).
    %
    % Output:
    %
    %   strategy - backtestStrategy object with properties that correspond
    %       to the parameters detailed above.  The backtestStrategy object
    %       is used in conjunction with the backtestEngine class to run
    %       backtests of portfolio investment strategies on historical
    %       data.
    %
    % Example:
    %
    %    % Load equity adjusted price data and convert to timetable
    %    T = readtable('dowPortfolio.xlsx');
    %    pricesTT = table2timetable(T(:,[1 3:end]),'RowTimes','Dates');
    %
    %    % Create backtest strategy
    %    rebalanceFcn = @(weights,priceWindow) ones(1,numel(weights)) / numel(weights);
    %    equalWeightStrategy = backtestStrategy("EqualWeight",rebalanceFcn,...
    %        'RebalanceFrequency',10,'TransactionCosts',[0.005 0.0025]);
    %
    %    % Create backtest engine
    %    backtester = backtestEngine(equalWeightStrategy);
    %
    %    % Run backtest and see results
    %    backtester = runBacktest(backtester,pricesTT);
    %    backtester.summary()
    %
    %   See also BACKTESTENGINE.
    
    % Copyright 2020 The MathWorks, Inc.
    properties
        Name
        RebalanceFcn
        RebalanceFrequency
        RebalanceDate
        TransactionCosts
        LookbackWindow
        InitialWeights
        
        % Fee input
        ManagementFeeRate
        IncentiveFeeRate
        %ManagementFeeFreq
        %IncentiveFeeFreq
        ManagementFeeDate
        IncentiveFeeDate
        HighWatermark
        HurdleFeeFreq
        Hurdle
        HurdleFeeRate

    end
    
    properties (SetAccess=protected, Hidden=true)
        Weights
    end
    
    methods
        function obj = backtestStrategy(name, rebalanceFcn, varargin)
            
            if nargin < 1
                name = "DefaultStrategy";
                rebalanceFcn = @(w,~,~) w;
            end
            
            obj.Name         = name;
            obj.RebalanceFcn = rebalanceFcn;
            
            ip = inputParser;
            ip.addParameter('RebalanceFrequency', 1);
            ip.addParameter('RebalanceDate', []);
            ip.addParameter('TransactionCosts', 0);
            ip.addParameter('LookbackWindow', [0 Inf]);
            ip.addParameter('InitialWeights', []);
            ip.addParameter('ManagementFeeRate', 0);
            ip.addParameter('IncentiveFeeRate', 0);
            %ip.addParameter('ManagementFeeFreq', 0);
            %ip.addParameter('IncentiveFeeFreq', 0);
            ip.addParameter('ManagementFeeDate', []);
            ip.addParameter('IncentiveFeeDate', []);

            ip.addParameter('HighWatermark', 0);
            ip.addParameter('HurdleFeeFreq', 0);
            ip.addParameter('Hurdle', 0);
            ip.addParameter('HurdleFeeRate', 0);
            
            
            
            ip.parse(varargin{:});
            result = ip.Results;
            
            obj.RebalanceFrequency = result.RebalanceFrequency;
            obj.RebalanceDate = result.RebalanceDate;
            obj.TransactionCosts   = result.TransactionCosts;
            obj.LookbackWindow     = result.LookbackWindow;
            obj.InitialWeights     = result.InitialWeights;
            % Fee
            obj.ManagementFeeRate = result.ManagementFeeRate;
            obj.IncentiveFeeRate = result.IncentiveFeeRate; 
            %obj.ManagementFeeFreq = result.ManagementFeeFreq;
            %obj.IncentiveFeeFreq = result.IncentiveFeeFreq;
            obj.ManagementFeeDate = result.ManagementFeeDate;
            obj.IncentiveFeeDate = result.IncentiveFeeDate;

            obj.HighWatermark = result.HighWatermark;
            obj.HurdleFeeFreq = result.HurdleFeeFreq;
            obj.Hurdle = result.Hurdle;
            obj.HurdleFeeRate = result.HurdleFeeRate;

        end
        
        % Property Setters
        function obj = set.Name(obj,value)
            % Must be a string or character vector
            validateattributes(value,["char","string"],"nonempty",mfilename,"Name");
            % Must be string and convertible to a valid variable name
            value = strrep(string(value)," ","_");
            if ~isvarname(value)
                error(message('finance:backtest:NoSpecialChars'));
            end
            obj.Name = value;
        end
        
        function obj = set.RebalanceFcn(obj,value)
            % Must be a function handle
            validateattributes(value,"function_handle","nonempty",mfilename,"RebalanceFcn");
            in  = nargin(value);
            out = nargout(value);
            validIn  = in < 0 || in == 2 || in == 3;
            validOut = out < 0 || out == 1;
            if validIn && validOut
                obj.RebalanceFcn = value;
            else
                error(message('finance:backtest:InvalidRebalanceFcn'));
            end
        end
        
        function obj = set.RebalanceFrequency(obj,value)
            % Must be a numeric scalar
            validateattributes(value,"numeric",["nonnegative" "scalar"],mfilename,"RebalanceFrequency");
            obj.RebalanceFrequency = value;
        end
        
        
        function obj = set.RebalanceDate(obj,value)
            % Must be a datetime vector
            
            if ~isempty(value)
                validateattributes(value,"datetime","vector",mfilename,"RebalanceDate");
            
                if ~all(isbusday(value))
                    error(message('finance:backtest:InvalidRebalanceDate Not Business Days'));
                else
                    obj.RebalanceDate = value;
                end
            else
                obj.RebalanceDate = value;
            end
        end
        
        function obj = set.ManagementFeeDate(obj,value)
            % Must be a datetime vector
            if ~isempty(value)
                validateattributes(value,"datetime","vector",mfilename,"ManagementFeeDate");
            
                if ~all(isbusday(value))
                    error(message('finance:backtest:InvalidManagementFeeDate Not Business Days'));
                else
                    obj.ManagementFeeDate = value;
                end
            else
                obj.ManagementFeeDate = value;
            end
        end
        
        function obj = set.IncentiveFeeDate(obj,value)
            % Must be a datetime vector
            if ~isempty(value)
                validateattributes(value,"datetime","vector",mfilename,"IncentiveFeeDate");
            
                if ~all(isbusday(value))
                    error(message('finance:backtest:InvalidIncentiveFeeDate Not Business Days'));
                else
                    obj.IncentiveFeeDate = value;
                end
            else
                obj.IncentiveFeeDate = value;
            end
        end

        function obj = set.ManagementFeeRate(obj,value)
            % Must be a numeric scalar
            validateattributes(value,"numeric","scalar",mfilename,"ManagementFeeRate");
            obj.ManagementFeeRate = value;
        end

        function obj = set.IncentiveFeeRate(obj,value)
            % Must be a numeric scalar
            validateattributes(value,"numeric","scalar",mfilename,"IncentiveFeeRate");
            obj.IncentiveFeeRate = value;
        end

        
        
        function obj = set.TransactionCosts(obj,value)
            % Must be a 1 or 2 element vector or else a function handle
            if isnumeric(value)
                if ismember(numel(value),[1 2])
                    obj.TransactionCosts = value(:)';
                else
                    error(message('finance:backtest:InvalidNumericTransactionCosts'));
                end
            elseif isa(value,'function_handle')
                obj.TransactionCosts = value;
            else
                error(message('finance:backtest:InvalidTransactionCosts'));
            end
        end
        
        function obj = set.LookbackWindow(obj,value)
            % Must be 1 or 2 element nonnegative sorted vector
            validateattributes(value,"numeric",["vector","nonnegative"],mfilename,"LookbackWindow");
            if ~ismember(numel(value),[1 2])
                error(message('finance:backtest:LookbackWindowSize'));
            end
            if sort(value(:)) ~= value(:)
                error(message('finance:backtest:LookbackWindowOrder'));
            end
            obj.LookbackWindow = value(:)';
        end
        
        function obj = set.InitialWeights(obj,value)
            % Must be empty or a numeric vector
            if isnumeric(value) && (isempty(value) || isvector(value))
                obj.InitialWeights = value(:)';
            else
                error(message('finance:backtest:InvalidInitialWeights'));
            end
        end
        
        function obj = set.Weights(obj,value)
            % Must be a numeric vector
            validateattributes(value,"numeric","vector",mfilename,"Weights");
            obj.Weights = value(:)';
        end
        
        function obj = rebalance(obj, currentWeights, prices, signal)
            % Rebalance portfolio weights.
            %
            % Syntax:
            %
            %   obj = rebalance(obj, currentWeights, prices)
            %   obj = rebalance(obj, currentWeights, prices, signal)
            %
            % Description:
            %
            %   The REBALANCE method is called by the backtestEngine to
            %   rebalance the portfolio positions of the backtestStrategy
            %   object using the function handle stored in the RebalanceFcn
            %   property.  See the help for the backtestStrategy class for
            %   details on how to specify the RebalanceFcn property.
            
            if nargin < 4
                obj.Weights = obj.RebalanceFcn(currentWeights, prices);
            else
                obj.Weights = obj.RebalanceFcn(currentWeights, prices, signal);
            end
        end
        
        function [buy, sell] = computeTransactionCosts(obj, deltaPositions)
            % Compute transaction costs from changes in asset positions.
            %
            % Syntax:
            %
            %   [buy, sell] = computeTransactionCosts(obj, deltaPositions)
            %
            % Description:
            %
            %   The COMPUTETRANSACTIONCOSTS method is called by the
            %   backtestEngine to compute the transaction costs incurred
            %   during a rebalance of the portfolio positions of the
            %   backtestStrategy object using the costs set in the
            %   TransactionCosts property.  See the help for the
            %   backtestStrategy class for details on how to specify the
            %   TransactionCosts property.
            
            if isa(obj.TransactionCosts,'function_handle')
                [buy, sell] = obj.TransactionCosts(deltaPositions);
            else
                rate = obj.TransactionCosts;
                [buy, sell] = percentTransactionCosts(deltaPositions, rate);
            end
            
        end
        
        function [AboveHurdleCost, UpdateLastPosition] = computePerformanceCost(obj, boolHurdle, currPortVal, updatelastposition)
            AboveHurdleCost = 0;
            currPortTotalVal = sum(currPortVal);
            UpdateLastPosition = updatelastposition;
            currReturn = (currPortTotalVal - UpdateLastPosition) /  UpdateLastPosition;
            if currReturn > obj.Hurdle
                AboveHurdleCost = (currReturn - obj.Hurdle) * UpdateLastPosition * obj.HurdleFeeRate;
                if boolHurdle
                    UpdateLastPosition = currPortTotalVal;
                end
            end
        end
        
        function [IncentiveCost,Update_HighWaterMark] = computeIncentiveCost(obj, highwatermark,bool_hwm, currPortVal,Mfee,NextIncFeeDate,PrevIncFeeDate)
            % highwatermark: a scalar that records the highest portfolio
            % value during the backtesting.
            %
            % currPortVal: a vector of positions (in market value), instead
            % of directly passing in the total mkt value of the portfolio.
            % Because a user may define their performance cost that is
            % asset specific. 
            
            IncentiveCost = 0;
            Update_HighWaterMark = highwatermark;
            currPortTotalVal = sum(currPortVal);
            if bool_hwm
                if  currPortTotalVal> highwatermark
                    IncentiveCost = (currPortTotalVal - highwatermark)*obj.IncentiveFeeRate * daysact(PrevIncFeeDate,  NextIncFeeDate)/360;
                    Update_HighWaterMark = currPortTotalVal; 
                end
            else 
                if currPortTotalVal> highwatermark+Mfee
                    IncentiveCost = (currPortTotalVal - highwatermark - Mfee)*obj.IncentiveFeeRate * daysact(PrevIncFeeDate,  NextIncFeeDate)/360;
                end
               
            end
            
        end
        
        function managementCost = computeManagementCost(obj, startPortfVal,NextMFeeDate,PrevMFeeDate)
            % endPortfVal: a vector of positions at the end of current trading 
            % (in market value), instead of directly passing in the total mkt value 
            % of the portfolio. Because a user may define their performance cost
            % that is asset specific. 
            managementCost = sum(startPortfVal)*obj.ManagementFeeRate*daysact(PrevMFeeDate,  NextMFeeDate)/360;
        end

    end
end

function [buy, sell] = percentTransactionCosts(deltaMktValue,rate)
% Flat percentile for transaction costs

if isscalar(rate)
    rate = [rate rate];
end
buy  = sum(max( deltaMktValue(:), 0)) * rate(1);
sell = sum(max(-deltaMktValue(:), 0)) * rate(2);

end
