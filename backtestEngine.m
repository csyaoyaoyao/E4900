classdef backtestEngine
    %BACKTESTENGINE Engine for portfolio backtesting.
    %
    % Syntax:
    %
    %   backtester = backtestEngine(strategies)
    %   backtester = backtestEngine(strategies, param1, value1,...)
    %
    % Description:
    %
    %   BACKTESTENGINE takes as input a vector of BACKTESTSTRATEGY objects
    %   and returns a BACKTESTENGINE object.  The backtest engine is used
    %   to backtest the portfolio trading strategies defined in the
    %   BACKTESTSTRATEGY objects.
    %
    % Input arguments:
    %
    %   strategies - A vector of backtestStrategy objects to be tested
    %     using the runBacktest() method.  Each backtestStrategy object
    %     defines a portfolio trading strategy.  See the help for
    %     backtestStrategy for more information.
    %
    % Optional Input Parameter Name/Value Pairs:
    %
    %   'RiskFreeRate' : A scalar that specifies the rate of return on
    %     uninvested capital in the portfolio (cash).  The RiskFreeRate is
    %     a decimal percentage and represents the risk free rate for one
    %     time step in the backtest.  For example, if the backtest is using
    %     daily asset price data, then the RiskFreeRate must be the daily
    %     rate of return for cash.  The default is 0.
    %
    %   'CashBorrowRate' : A scalar that specifies the rate of interest
    %     accrual on negative cash balances (margin) during the backtest.
    %     The CashBorrowRate is a decimal percentage and represents the
    %     interest accrual rate for one time step in the backtest.  For
    %     example, if the backtest is using daily asset price data, then
    %     the CashBorrowRate must be the daily interest rate for negative
    %     cash balances.  The default is 0.
    %
    %   'InitialPortfolioValue' : A scalar that specifies the initial
    %     portfolio value for all strategies.  The default is 10,000.
    %
    % Output:
    %
    %   backtester - backtestEngine object with properties that correspond
    %       to the parameters detailed above.  The backtestEngine object is
    %       used to run backtests of portfolio investment strategies on
    %       historical data.
    %
    %   Additionally the backtestEngine object has the following properties
    %   which are empty until the backtest is run using the runBacktest
    %   method:
    %
    %     * NumAssets : The number of assets in the portfolio universe.
    %       Derived from the timetable of adjusted prices passed into the
    %       runBacktest method.
    %
    %     * Returns : A NumTimeSteps-by-NumStrategies timetable of strategy
    %       returns.  Returns are per time step.  For example if daily
    %       prices are used in the runBacktest method then Returns will be
    %       the daily strategy returns.
    %
    %     * Positions : A struct containing a NumTimeSteps-by-NumAssets
    %       timetable of asset positions for each strategy.  For example if
    %       daily prices are used in the runBacktest method then the
    %       Positions struct will hold timetables containing the daily
    %       asset positions.
    %
    %     * Turnover : A NumTimeSteps-by-NumStrategies timetable of
    %       strategy turnover.
    %
    %     * BuyCost : A NumTimeSteps-by-NumStrategies timetable of
    %       transaction costs for the asset purchases of each strategy.
    %
    %     * SellCost : A NumTimeSteps-by-NumStrategies timetable of
    %       transaction costs for the asset sales of each strategy.
    %
    % Example:
    %
    %    % Load equity adjusted price data and convert to timetable
    %    T = readtable('dowPortfolio.xlsx');
    %    pricesTT = table2timetable(T(:,[1 3:end]),'RowTimes','Dates');
    %
    %    % Create backtest strategy
    %    rebalanceFcn = @(weights,priceWindow) ones(1,numel(weights)) / numel(weights);
    %    equalWeightStrategy = backtestStrategy("EqualWeight",rebalanceFcn);
    %
    %    % Create backtest engine
    %    backtester = backtestEngine(equalWeightStrategy,...
    %        'RiskFreeRate',0.2,...
    %        'InitialPortfolioValue',1e6);
    %
    %    % Run backtest and see results
    %    backtester = runBacktest(backtester,pricesTT);
    %    backtester.summary()
    %
    %   See also BACKTESTSTRATEGY.
    
    % Copyright 2020 The MathWorks, Inc.
    properties
        Strategies
        RiskFreeRate
        CashBorrowRate
        InitialPortfolioValue
    end
    
    properties (SetAccess=private)
        NumAssets
        Returns
        Positions
        TWReturn
        MDReturn
        Turnover
        BuyCost
        SellCost
        PwgtStructEOD
        Fee
    end
    
    methods
        function obj = backtestEngine(strategies, varargin)
            
            if nargin < 1
                strategies = backtestStrategy;
            end
            
            obj.Strategies = strategies;
            
            ip = inputParser;
            ip.addParameter('RiskFreeRate', 0);
            ip.addParameter('CashBorrowRate', 0);
            ip.addParameter('InitialPortfolioValue', 1e4);
            
            ip.parse(varargin{:});
            result = ip.Results;
            
            obj.RiskFreeRate          = result.RiskFreeRate;
            obj.CashBorrowRate        = result.CashBorrowRate;
            obj.InitialPortfolioValue = result.InitialPortfolioValue;
        end
        
        % Property Setters
        function obj = set.Strategies(obj,value)
            % Must be a vector of backtestStrategy objects
            validateattributes(value,"backtestStrategy","vector",mfilename,"Strategies");
            % All non-empty InitialWeights vectors must have same size
            initialWeights = {value.InitialWeights};
            nonEmptyIdx = ~cellfun(@isempty,initialWeights);
            initialWeightsSz = cellfun(@size,initialWeights,'UniformOutput',false);
            specifiedSizes = initialWeightsSz(nonEmptyIdx);
            if ~isempty(specifiedSizes)
                size1 = specifiedSizes{1};
                if ~all(cellfun(@(si) isequal(si, size1), specifiedSizes))
                    error(message('finance:backtest:StrategySizeMismatch'));
                end
            end
            % All strategies must have unique names
            names = [value.Name];
            if numel(names) ~= numel(unique(names))
                error(message('finance:backtest:NonUniqueStrategies'));
            end
            obj.Strategies = value(:)';
        end
        
        
        function obj = set.RiskFreeRate(obj,value)
            % Must be a numeric scalar
            validateattributes(value,"numeric","scalar",mfilename,"RiskFreeRate");
            obj.RiskFreeRate = value;
        end

        
        function obj = set.CashBorrowRate(obj,value)
            % Must be a numeric scalar
            validateattributes(value,"numeric","scalar",mfilename,"CashBorrowRate");
            obj.CashBorrowRate = value;
        end
        
        function obj = set.InitialPortfolioValue(obj,value)
            % Must be a positive numeric scalar
            validateattributes(value,"numeric",["scalar","positive"],mfilename,"InitialPortfolioValue");
            obj.InitialPortfolioValue = value;
        end
        
        function obj = runBacktest(obj, pricesTT, varargin)
            % Run backtest.
            %
            % Syntax:
            %
            %   backtester = runBacktest(backtester, pricesTT)
            %   backtester = runBacktest(backtester, pricesTT, param1, value1,...)
            %
            %   backtester = runBacktest(backtester, pricesTT, signalTT)
            %   backtester = runBacktest(backtester, pricesTT, signalTT, param1, value1,...)
            %
            % Description:
            %
            %   The runBacktest method runs the backtest of the strategies
            %   using the adjusted asset price data as well as the signal
            %   data if provided.  After completion, runBacktest will
            %   populate several properties in the backtestEngine object
            %   with the results of the backtest.
            %
            % Input arguments:
            %
            %   pricesTT - Asset prices.  A timetable of asset prices used
            %     by the backtestEngine to backtest the strategies.  Must
            %     be specified as a timetable where each column contains a
            %     time series of prices for an investable asset. Historical
            %     prices should be adjusted for splits and dividends.
            %
            %   signalTT - Signal data.  A timetable of trading signal data
            %     that the strategies can use to make trading decisions.
            %     The signal argument is optional.  If provided, the
            %     backtestEngine will call the strategy rebalance functions
            %     with both asset price data as well as signal data.  The
            %     signal data timetable must have the same time dimension
            %     as the asset price timetable.
            %
            %   Transaction - Transaction data. A timetable of transactions,
            %     The transaction argument is optional.  If provided, the
            %     backtestEngine will extract the daily cashflow from the
            %     data to help calculate time-weighted and MD returns. The
            %     transaction data timetable must have each date for each row
            %     and each asset in the portfolio for each column.
            %
            % Optional Input Parameter Name/Value Pairs:
            %
            %   'Start' : A scalar that sets the starting time for the
            %     backtest.  The Start parameter is specified in terms of
            %     time steps, i.e. the row of the asset price timetable
            %     where the backtest is to begin.
            %
            %   'End' : A scalar that sets the ending time for the
            %     backtest.  The End parameter is specified in terms of
            %     time steps, i.e. the final row of the asset price
            %     timetable that is included in the backtest.
            %
            %  The runBacktest method initializes each strategy to the
            %  InitialPortfolioValue and begins walking through the
            %  asset price data.
            %
            %  At each time step the engine applies the asset returns to
            %  all assets positions. The engine then determines if the
            %  strategy needs to be rebalanced based on the
            %  RebalanceFrequency and LookbackWindow properties of the
            %  backtestStrategy objects.
            %
            %  For strategies that need rebalancing, the backtestEngine
            %  calls their rebalance functions with a rolling window of
            %  asset price data (and signal data if provided) based on the
            %  LookbackWindow property of each strategy.
            %
            %  Transaction costs are then paid based on the changes in
            %  asset positions and the TransactionCosts property of the
            %  backtestStrategy objects.
            %
            %  After completing the backtest, the results are stored in
            %  several properties.  See the help for backtestEngine and
            %  backtestStrategy for more details.
            
            % Verify we have a valid prices timetable
            if ~isa(pricesTT,'timetable')
                error(message('finance:backtest:InvalidAssetTT'));
            end
            
            % Parse out and validate Start and End times before checking
            % for missing data.
            ip = inputParser;
            ip.addOptional('SignalTT', []);
            ip.addOptional('Transaction', []);

            ip.addParameter('Start', 1,...
                @(p) validateattributes(p,"numeric","scalar",...
                "runBacktest","Start"));
            ip.addParameter('End', size(pricesTT,1),...
                @(p) validateattributes(p,"numeric","scalar",...
                "runBacktest","End"));
            
            ip.parse(varargin{:});
            result = ip.Results;
            
            signalTT  = result.SignalTT;
            %Transaction  = result.Transaction;
            startTime = result.Start;
            endTime   = result.End;
            
            % Validate Start & End
            if startTime < 1 || endTime <= startTime || size(pricesTT,1) <= startTime
                error(message('finance:backtest:InvalidStart'));
            end
            if endTime <= 1 || size(pricesTT,1) < endTime
                error(message('finance:backtest:InvalidEnd'));
            end
            
            % Validate asset prices in backtest range
            missingAssetPrices = ismissing(pricesTT(startTime:endTime,:));
            if any(missingAssetPrices(:))
                error(message('finance:backtest:MissingAssetData'));
            end
            assetTimeDim = pricesTT.Properties.DimensionNames{1};
            
            % Validate signal data
            if isempty(signalTT)
                signalSpecified = false;
            else
                signalSpecified = true;
                if ~isa(signalTT,'timetable')
                    error(message('finance:backtest:InvalidSignalTT'));
                end
                % Verify sizes match before checking for missing values
                % since start and end times were validated against the
                % prices timetable
                if size(pricesTT,1) ~= size(signalTT,1)
                    error(message('finance:backtest:SignalSizeMismatch'));
                end
                missingSignals = ismissing(signalTT(startTime:endTime,:));
                if any(missingSignals(:))
                    error(message('finance:backtest:MissingSignalData'));
                end
                assetTime = pricesTT.(assetTimeDim);
                signalTimeDim = signalTT.Properties.DimensionNames{1};
                signalTime = signalTT.(signalTimeDim);
                if ~isequal(assetTime,signalTime)
                    error(message('finance:backtest:SignalTimeMismatch'));
                end
            end
            
            % backtest time steps are in terms of returns, not prices
            backtestRange = startTime:endTime-1;
            numTimeSteps  = numel(backtestRange);
            
            assetReturns  = tick2ret(pricesTT);  % in percentage
            obj.NumAssets = size(assetReturns.Variables,2)-1;
            assetIdx      = find(~strcmp("overdraft",assetReturns.Properties.VariableNames)); %1:obj.NumAssets;
            
            strategies  = obj.Strategies;
            nStrategies = numel(strategies);
            
            % Variables to store intermediary results
            returns  = zeros(numTimeSteps, nStrategies);
            turnover = zeros(numTimeSteps, nStrategies);
            buycost  = zeros(numTimeSteps, nStrategies);
            sellcost = zeros(numTimeSteps, nStrategies);
            
            initialPositions = zeros(nStrategies, obj.NumAssets);
            
            % Run backtest for each strategy
            for i = 1:numel(strategies)
                
                % Initialize position array
                stratName = strategies(i).Name;
                positionsEOD.(stratName) = zeros(numTimeSteps, obj.NumAssets );
                
                TW_Return_Full.(stratName) = zeros(numTimeSteps, obj.NumAssets);  % all assets
                MD_Return_Full.(stratName) = zeros(numTimeSteps, obj.NumAssets);  % all assets
                eodWeights_Full.(stratName) = zeros(numTimeSteps, obj.NumAssets);
                fee.(stratName) = zeros(numTimeSteps, 4);
                
                % Initialize highwatermark to InitialPortfolioValue
                highwatermark = obj.InitialPortfolioValue;

                updatelastposition = obj.InitialPortfolioValue;

                % Set initial positions
             
                if ~isempty(strategies(i).InitialWeights)
                    
                    % If the strategy specifies initial weights, set them
                    if numel(strategies(i).InitialWeights) ~= obj.NumAssets
                        error(message('finance:backtest:InvalidInitialWeightsSize'));
                    end
                    initialPositions(i,:) = obj.InitialPortfolioValue * [strategies(i).InitialWeights];
                else
                    % Otherwise, if no weights are set, we set it equal to first line of signalTT
                    strategies(i).InitialWeights = signalTT{1,assetIdx}/sum(signalTT{1,assetIdx});
                    initialPositions(i,:) = obj.InitialPortfolioValue * [strategies(i).InitialWeights];
                end
                
                % Run backtest
                cash_idx = find(pricesTT.Properties.VariableNames=="cash");
                overdraft_idx = find(pricesTT.Properties.VariableNames=="overdraft");
                
                % Initiate the inputs for Modified Dietz Method and Time
                % weighted Return 
                deltaCashFlow_array = []; %zeros(1,obj.NumAssets+1);
                Position_start = initialPositions(i,(1:end));
                prevPositions = initialPositions(i,:);
                
                numRebalance = length(strategies(i).RebalanceDate);
                
                % Initialize Rebalance Dates
                if ~isempty(strategies(i).RebalanceDate)
                    strategies(i).RebalanceDate = unique(strategies(i).RebalanceDate);
                    NextRebalanceDate = strategies(i).RebalanceDate(1);
                    Rebalance_anchor = 1;
                else
                    NextRebalanceDate = datestr(1);
                    Rebalance_anchor = 0;
                end
                
                %%% Initiate Fees
                % Initialize ManagementFee
                MFeePaid = 0;                
                numMFee = length(strategies(i).ManagementFeeDate);
                if ~isempty(strategies(i).ManagementFeeDate)
                    strategies(i).ManagementFeeDate = unique(strategies(i).ManagementFeeDate);
                    NextMFeeDate = strategies(i).ManagementFeeDate(1);
                    MFee_anchor = 1;
                else
                    NextMFeeDate = datestr(1);
                    MFee_anchor = 0;
                end
                
                % Initialize IncentiveFee
                numIncFee = length(strategies(i).IncentiveFeeDate);
                PrevIncFeeDate = pricesTT.Time(1);
                if ~isempty(strategies(i).IncentiveFeeDate)
                    strategies(i).IncentiveFeeDate = unique(strategies(i).IncentiveFeeDate);
                    NextIncFeeDate = strategies(i).IncentiveFeeDate(1);
                    IncFee_anchor = 1;
                else
                    NextIncFeeDate = datestr(1);
                    IncFee_anchor = 0;
                end
                
                % Boolean highwatermark
                Bool_HighWatermark = strategies(i).HighWatermark;
                
                
                for btIdx = 1:numTimeSteps
                    % time step returns
                    returnIdx     = backtestRange(btIdx);
                    assetReturnsi = assetReturns{returnIdx, assetIdx};
                    delta_position = zeros(1,obj.NumAssets);
                    
                    % end of time step price/signal index
                    endTimeIdx = backtestRange(btIdx) + 1;
                    
                    % Start of day positions
                    sodPositions = prevPositions;
                    sodPortValue = sum(sodPositions);
                    
                    % Daily returns
                    if 0 <= sodPositions(cash_idx)
                        % Zero or positive cash
                        retn = 1 + assetReturnsi;
                    else
                        % Negative cash
                        assetReturnsi(cash_idx) = assetReturns{returnIdx,overdraft_idx};
                        retn = 1 + assetReturnsi;
                    end
                    
                    % Apply current day's return
                    eodPositions = sodPositions .* retn;
                    eodPortValue = sum(eodPositions);
                    eodAssetWeights = eodPositions(assetIdx) / eodPortValue;
                    
                    
                    % Rebalance if necessary
                    if (pricesTT.Time(btIdx) ==  NextRebalanceDate)
                        if Rebalance_anchor<numRebalance
                            Rebalance_anchor = Rebalance_anchor+1;
                            NextRebalanceDate = strategies(i).RebalanceDate(Rebalance_anchor);
                        end
                        
                        % Verify we have enough data for lookback window
                        %if 1 <= endTimeIdx - strategies(i).LookbackWindow(1)
                        if 0 <= endTimeIdx - strategies(i).LookbackWindow(1)
                            % Generate rolling windows
                            startTimeIdx = endTimeIdx - strategies(i).LookbackWindow(end); %；-1  ; % + 1;
                            dataRange = max(1,startTimeIdx):(endTimeIdx-1);%-2;  %min(endTimeIdx,numTimeSteps);
                            assetLookbackData = pricesTT(dataRange, assetIdx); %assetReturns(dataRange, assetIdx); %%% pricesTT(dataRange, :);
                            
                            % Rebalance
                            if signalSpecified && nargin(strategies(i).RebalanceFcn) ~= 2
                                % Signal data was specified and rebalance
                                % function will take 3 inputs or varargin.
                                signalLookbackData = signalTT(dataRange, assetIdx);
                                %strategies(i) = rebalance(strategies(i), eodAssetWeights, assetLookbackData, eodAssetWeights);
                                strategies(i) = rebalance(strategies(i), eodAssetWeights, assetLookbackData, signalLookbackData);
                                
                            else
                                strategies(i) = rebalance(strategies(i), eodAssetWeights, assetLookbackData);
                            end
                            
                            eodAssetWeightsUpdated = strategies(i).Weights/ sum(strategies(i).Weights);
                            deltaWeights = eodAssetWeightsUpdated - eodAssetWeights;
                            
                            %
                            % Record daily metrics
                            turnover(btIdx,i) = sum(abs(deltaWeights)) / 2;
                            [buycost(btIdx,i), sellcost(btIdx,i)] = strategies(i).computeTransactionCosts(...
                                eodPortValue * deltaWeights(1:end ~= cash_idx));
                            
                            delta_cash = delta_cash_flow_vec(eodPortValue * deltaWeights,buycost(btIdx,i),sellcost(btIdx,i),cash_idx);
                            delta_position = eodPortValue * deltaWeights; 
                            eodPositionsUpdated = eodPositions + delta_cash; 
                            
                            eodPositionsUpdated(cash_idx) = eodPositionsUpdated(cash_idx) + sum(eodPositionsUpdated((1:end~=cash_idx)&(eodPositionsUpdated<0)));
                            eodPositionsUpdated((1:end~=cash_idx)&(eodPositionsUpdated<0))=0;
                            
                            %eodPortValueUpdated = eodPortValue - buycost(btIdx, i) - sellcost(btIdx, i);
                            
                            % Delta Cash Flow
                            deltaCashFlow = eodPositionsUpdated- eodPositions;

                            % Transaction fee is stored in second column
                            fee.(strategies(i).Name)(btIdx,2) = buycost(btIdx, i) + sellcost(btIdx, i);
                            
                            % Pay transaction fees & update weights
                            %eodPortValue = eodPortValueUpdated;
                            eodPortValue = sum(eodPositionsUpdated);
                            eodAssetWeights = eodPositionsUpdated/sum(eodPositionsUpdated);
                            eodPositions = eodPositionsUpdated;
                            % Positive weights should always be long positions
                            if eodPortValue < 0
                                eodAssetWeights = -eodAssetWeights;
                            end
                        end
                    end
                    
                    %if btIdx>=62
                    %    disp("ok");
                    %end
                    % Update position if any transaction took palce
                    % If transaction didn't take place
                    % Transaction_value(btIdx) is a zero vector, and it
                    % won't effect Position, value or weight
                    
                    %eodPositions = eodPositions + Transaction_value(btIdx,:);
                    %eodPortValue = eodPortValue + sum(Transaction_value(btIdx,:));
                    %eodPortValue = sum(eodPositions);
                    
                    
                    %%%Calculate Management & Incentive Fee
                    MFee = 0;
                    IncentiveFee = 0;
                    AboveHurdleCost = 0;
                    
                    % Calcualte Management Fee
                    if pricesTT.Time(btIdx) ==  NextMFeeDate
                        %MFee = strategies(i).computeManagementCost(StartPortfolioValue_withupdate,NextMFeeDate,PrevMFeeDate);
                        PrevMFeeDate = NextMFeeDate;
                        if MFee_anchor<numMFee
                            MFee_anchor = MFee_anchor+1;
                            NextMFeeDate = strategies(i).ManagementFeeDate(MFee_anchor);
                        end
                       
                        MFee = strategies(i).computeManagementCost(sodPositions,NextMFeeDate,PrevMFeeDate);
                        
                        MFeePaid = MFee;
                        % StartPortfolioValue_withupdate = eodPortValue;
                        % Update HighWaterMark
                        if Bool_HighWatermark
                            highwatermark = eodPortValue;
                        end
                        
                    end
                    
                    % Calculate Incentive Fee
                    if pricesTT.Time(btIdx) ==  NextIncFeeDate 
                        eodMktVal = eodPortValue;
                        [IncentiveFee, HighWaterMark]= strategies(i).computeIncentiveCost(highwatermark,Bool_HighWatermark,eodMktVal,MFeePaid,NextIncFeeDate,PrevIncFeeDate);
                        highwatermark = HighWaterMark;
                        PrevIncFeeDate = NextIncFeeDate;
                        if IncFee_anchor<numIncFee
                            IncFee_anchor = IncFee_anchor+1;
                            NextIncFeeDate = strategies(i).IncentiveFeeDate(IncFee_anchor);
                        end                  
                    end

                    % Calculate Hurdle Rate
                    if mod(btIdx, strategies(i).HurdleFeeFreq) == 0
                        eodMktVal = eodPortValue;
                        boolHurdle = strategies(i).Hurdle;
                        [AboveHurdleCost, UpdateLastPosition] = strategies(i).computePerformanceCost(boolHurdle, eodMktVal, updatelastposition);
                        updatelastposition = UpdateLastPosition;
                    end
                    
                    fee.(strategies(i).Name)(btIdx,4) = IncentiveFee + MFee + AboveHurdleCost;
                    cash_position = eodPositions(cash_idx)- fee.(strategies(i).Name)(btIdx,2) - fee.(strategies(i).Name)(btIdx,4);
                    fee.(strategies(i).Name)(btIdx,1) = cash_position;
                    
                    if cash_position<0
                        fee.(strategies(i).Name)(btIdx,3) = -cash_position * assetReturns{btIdx,overdraft_idx};
                    end
                    
                    deltaCashFlow(cash_idx) = deltaCashFlow(cash_idx) - fee.(strategies(i).Name)(btIdx,2) - fee.(strategies(i).Name)(btIdx,4);
                    
                    % Record positions and returns after Paying all fees
                    eodPositions(cash_idx) = cash_position; 
                    prevPositions = eodPositions;

                    eodPortValue = sum(eodPositions);
                    
                    %%%%% 
                    %returns(btIdx, i) = eodPortValue / sodPortValue - 1;
                    %returns(btIdx, i) = (eodPortValue - sodPortValue-sum(deltaCashFlow))/ (sodPortValue-sum(deltaCashFlow));
                    returns(btIdx, i) = (eodPortValue - sodPortValue-sum(deltaCashFlow))/ (sodPortValue);
                    returns(btIdx, i) = (eodPortValue - sodPortValue)/ (sodPortValue);

                    positionsEOD.(strategies(i).Name)(btIdx, :) = eodPositions;
                    
                    % Record EOD weights
                    eodAssetWeights = eodPositions / eodPortValue;
                    eodWeights_Full.(strategies(i).Name)(btIdx, :) = eodAssetWeights;
                    
                    
                    % Calculate Time Weighted Return 
                    asset_ind_TW = find(sodPositions ~= 0); % If sodPositions==0, the return during this day is just 0.
                    TW_PNL = zeros(1,obj.NumAssets);
                    %TW_PNL(asset_ind_TW) = eodPositions(asset_ind_TW) - sodPositions(asset_ind_TW) - deltaCashFlow(asset_ind_TW);
                    TW_PNL(asset_ind_TW) = eodPositions(asset_ind_TW) - sodPositions(asset_ind_TW) - delta_position(asset_ind_TW);
                    TW_return = TW_PNL./sodPositions;
                    
                    % Record Time Weighted Return
                    TW_Return_Full.(strategies(i).Name)(btIdx, :) = TW_return;

                end
            end

            VarNames = [assetReturns.Properties.VariableNames(assetIdx)];
            
            % Positions uses the prices time since it includes initial T=0
            %Time = assetReturns.(assetTimeDim)(startTime:endTime);
            %assetReturns
            for i = 1:numel(strategies)
                stratName = strategies(i).Name;
                Time = pricesTT.(assetTimeDim)(startTime:endTime);

                % Add the initial positions to the position data
                positionData = [initialPositions(i,assetIdx); positionsEOD.(stratName)];
                positionTable = array2table(positionData, 'VariableNames', VarNames);
                positionTable = addvars(positionTable, Time, 'Before', VarNames(1));
                obj.Positions.(stratName) = table2timetable(positionTable);
                
                Time = assetReturns.(assetTimeDim)(startTime:endTime-1);
                % Return Time Weighted Return
                TWReturnData = [TW_Return_Full.(stratName)];
                TWReturnTable = array2table(TWReturnData, 'VariableNames', VarNames);
                TWReturnTable = addvars(TWReturnTable, Time,'Before',VarNames(1));
                obj.TWReturn.(stratName) = table2timetable(TWReturnTable);

                % Return Modified Dietz Return
                MDReturnData = [MD_Return_Full.(stratName)]; %zeros(1,obj.NumAssets);
                MDReturnTable = array2table(MDReturnData, 'VariableNames', VarNames);
                MDReturnTable = addvars(MDReturnTable, Time,'Before',VarNames(1));
                obj.MDReturn.(stratName) = table2timetable(MDReturnTable);
                
                % EOD Stucture Weighted                zeros(1,obj.NumAssets)
                eodWeightsData = [eodWeights_Full.(stratName)];
                eodWeightsTable = array2table(eodWeightsData, 'VariableNames', VarNames);
                eodWeightsTable = addvars(eodWeightsTable, Time,'Before',VarNames(1));
                obj.PwgtStructEOD.(stratName) = table2timetable(eodWeightsTable);

                feeData = [fee.(stratName)];
                feeTable = array2table(feeData, 'VariableNames', ["CashAccount", "BrokerFee", "BankFee", "ManagementFeeAccount"]);
                feeTable = addvars(feeTable, Time, 'Before', 'CashAccount');
                obj.Fee.(stratName) = table2timetable(feeTable);

            end
            
            % Remaining metrics use returns times
            Time = assetReturns.(assetTimeDim)(startTime:endTime-1);
            
            % Extract strategy names
            stratsName = [strategies.Name];
            
            returns = array2table(returns,'VariableNames',stratsName);
            returns = addvars(returns, Time,'Before',stratsName{1});
            obj.Returns = table2timetable(returns);
            
            turnover = array2table(turnover,'VariableNames',stratsName);
            turnover = addvars(turnover, Time,'Before',stratsName{1});
            obj.Turnover = table2timetable(turnover);
            
            buycost = array2table(buycost,'VariableNames',stratsName);
            buycost = addvars(buycost, Time,'Before',stratsName{1});
            obj.BuyCost = table2timetable(buycost);
            
            sellcost = array2table(sellcost,'VariableNames',stratsName);
            sellcost = addvars(sellcost, Time,'Before',stratsName{1});
            obj.SellCost = table2timetable(sellcost);
            

        end
        
        
        function summaryTable = summary(obj)
            % Generate summary table of backtest results.
            %
            % Syntax:
            %
            %   summaryTable = summary(backtester)
            %
            % Description:
            %
            %   The summary method generates a table of metrics to
            %   summarize the backtest.  Each row of the table is a
            %   calculated metric and each column represents a strategy.
            %
            %   The reported metrics are:
            %
            %   TotalReturn - The total return of the strategy over the
            %       entire backtest.
            %   SharpeRatio - The Sharpe ratio for each strategy.
            %   Volatility - The volatility of each strategy over the
            %       backtest.
            %   AverageTurnover - Average turnover per time step as a
            %       decimal percent.
            %   MaxTurnover - Maximum turnover in a single time step.
            %   AverageReturn - Average return per time step.
            %   MaxDrawdown - Maximum portfolio drawdown as a decimal percent.
            %   AverageBuyCost - Average per time step transaction costs
            %       for asset purchases.
            %   AverageSellCost - Average per time step transaction costs
            %       for asset sales.
            
            if isempty(obj.Positions)
                error(message('finance:backtest:NoData'));
            end
            
            % Compute metrics summary table
            dailyReturns    = obj.Returns.Variables;
            averageReturns  = mean(dailyReturns);
            stdReturns      = std(dailyReturns);
            sharpeRatio     = sharpe(dailyReturns,obj.RiskFreeRate/360);
            compoundReturns = ret2tick(dailyReturns, 'method', 'simple');
            
            % maxdrawdown does not support non-positive portfolio values
            negativeIdx = any(compoundReturns <= 0);
            if any(negativeIdx)
                maxDrawdown = nan(size(sharpeRatio));
                maxDrawdown(~negativeIdx) = maxdrawdown(compoundReturns(:,~negativeIdx));
            else
                maxDrawdown = maxdrawdown(compoundReturns);
            end
            
            totalReturn     = compoundReturns(end, :)-1;
            avgTurnover     = mean(obj.Turnover.Variables);
            avgBuyCost      = mean(obj.BuyCost.Variables);
            avgSellCost     = mean(obj.SellCost.Variables);
            maxTurnover     = max(obj.Turnover.Variables);
            
            metricData = [totalReturn', sharpeRatio', stdReturns',...
                avgTurnover', maxTurnover', averageReturns', ...
                maxDrawdown', avgBuyCost', avgSellCost']';
            
            metricNames = ["TotalReturn",  "SharpeRatio", "Volatility",...
                "AverageTurnover", "MaxTurnover", "AverageReturn", ...
                "MaxDrawdown", "AverageBuyCost", "AverageSellCost"];
            
            colStratsNames = (obj.Returns.Properties.VariableNames)';
            
            summaryTable = array2table(metricData, 'VariableNames', colStratsNames, 'RowNames', metricNames);
            
        end
    end
end


 
function delta_cash_flow = delta_cash_flow_vec(deltaPositions,buy,sell,cash_idx)
% find delta cash flow
delta_cash_flow = zeros(1,length(deltaPositions));
pos_idx = find(deltaPositions>0);
neg_idx = find(deltaPositions<0);
neg_idx = neg_idx(neg_idx~=cash_idx);
delta_cash_flow(pos_idx) = deltaPositions(pos_idx)- buy*deltaPositions(pos_idx)/sum(deltaPositions(pos_idx));
delta_cash_flow(neg_idx) = deltaPositions(neg_idx)- sell*deltaPositions(neg_idx)/sum(deltaPositions(neg_idx));
delta_cash_flow(cash_idx) = deltaPositions(cash_idx);

end
