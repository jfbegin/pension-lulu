classdef DCALMStudy < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ProjYear
        nscen
        z0
        InfAnnual
        SrAnnual
        DivAnnual
        BondAnnual
        BondYield15
        BondYield1m
        
        EntryAge
        RetAge
        radix
        w
        nStart
        
        RetYears
        rates
        AnnFactors
        
        StartSal
        meric
        HistInf
        
        salary0
        StartPayroll
        CumInf
        salary
        
%         StockWeight
%         BondWeight
%         ValRateTypeContri
%         EquityRPContri
%         LboundContri
%         UboundContri
%         FixValRateContri
        
        TargetBen
        ContriRate
        
%         ValRateType
%         EquityRP
%         Lbound
%         Ubound
%         FixValRate
        ValRateMat
        InitialValRate
        
        ARMat
%        InitialLiaPUC %InitialLiaPUC
        InitialFund
        InitialFundGen
        
%        InitialFundRatio
%         NoActionRange
%         TriggerLFR
%         TriggerUFR
%         AdjAction
%         AdjLFR
%         AdjUFR
        
        SimulationResult
        
        mt
        Mt
        
        % struct
        TargetCalAssump
        DecumAssump
        PlanDesign
        
        ReplaceRatioAtRetire
        
        DistReturnList
    end
    
    methods
        function obj = DCALMStudy(ProjYear, nscen, assetScenario, TargetCalAssump, PlanDesign)
           %% Read in AssetInput
            obj.ProjYear = ProjYear; %55
            obj.nscen = nscen; %7788
            %continuously compounded, get from the last column of the hist
            %data
%            obj.z0 = [-7.7464; -6.4462; 0.0023; 0.0382676; -5.9976];  
            obj.z0 = [-6.34;-5.6802;0.0017; 0.0008; -6.3029];
            obj.z0([1 2 5]) = exp(obj.z0([1 2 5]));
            obj.PlanDesign = PlanDesign;
%             obj.InfAnnual = csvread("../../AssetInput/infAnnualEfft.csv");
%             obj.SrAnnual = csvread("../../AssetInput/srAnnualEfft.csv");
%             obj.DivAnnual = csvread("../../AssetInput/divAnnualEfft.csv");
%             obj.BondAnnual = csvread("../../AssetInput/bondAnnualEfft.csv");
%             obj.BondYield15 = csvread("../../AssetInput/bondYield_L_Efft.csv");
%             obj.BondYield1m = csvread("../../AssetInput/bondYield_S_Efft.csv");

%             obj.InfAnnual = round(csvread("../../Data/AssetInput/infAnnualEfft.csv"),6);
%             obj.SrAnnual = round(csvread("../../Data/AssetInput/srAnnualEfft.csv"),6);
%             obj.DivAnnual = round(csvread("../../Data/AssetInput/divAnnualEfft.csv"),6);
%             obj.BondAnnual = round(csvread("../../Data/AssetInput/bondAnnualEfft.csv"),6);
%             obj.BondYield15 = round(csvread("../../Data/AssetInput/bondYield_L_Efft.csv"),6);
%             obj.BondYield1m = round(csvread("../../Data/AssetInput/bondYield_S_Efft.csv"),6);

            obj.InfAnnual = assetScenario.infAnnualEfft;
            obj.SrAnnual = assetScenario.srAnnualEfft;
            obj.DivAnnual = assetScenario.divAnnualEfft;
            obj.BondAnnual = assetScenario.bondAnnualEfft;
            obj.BondYield15 = assetScenario.bondYield_L_Efft;
            obj.BondYield1m = assetScenario.bondYield_S_Efft;
            
            obj.mt = assetScenario.mt;
            obj.Mt = assetScenario.Mt;
           %% Set Up - demographic variables
            obj.EntryAge = 30;
            obj.RetAge = 65;
            obj.radix = 100;
            
           %% Simplified Mortality Table
            obj.w = 85;
            obj.nStart = zeros(obj.w, 1);
            obj.nStart(obj.EntryAge:obj.w) = ones((obj.w - obj.EntryAge + 1), 1) .* 100;
            
           %% Annuity Factor with fix valuation Rate
            %Calculate annuity factors recursively for each possible valuation rate from 1% to 25%
            %Assume annual payment at BOY            
            obj.RetYears = obj.w - obj.RetAge;
            obj.rates = (0.01:0.001:0.250)';
            obj.AnnFactors = zeros(obj.w, size(obj.rates, 1));
            obj.AnnFactors(obj.w, :) = 1;
            for i=(obj.w-1):-1:obj.RetAge 
                obj.AnnFactors(i,:) = obj.AnnFactors(i + 1,:) .* obj.nStart(i + 1) ./ obj.nStart(i) ./ (1 + obj.rates') + 1;
            end
            for i=(obj.RetAge-1):-1:obj.EntryAge 
                obj.AnnFactors(i,:) = obj.AnnFactors(i + 1,:) .* obj.nStart(i + 1) ./ obj.nStart(i) ./ (1 + obj.rates');
            end
            obj.AnnFactors(1:(obj.EntryAge - 1),:) = 0;
            
%             obj.ContriAnnFactors = zeros((obj.w-obj.EntryAge +1), size(obj.rates,1));
% 
%             for i = 1:1:(obj.RetAge - obj.EntryAge)
%                 obj.ContriAnnFactors(i,:) = (1-(1./(1+obj.rates').^(obj.RetAge - obj.EntryAge-i+1)))./(1-(1./(1+obj.rates')));
%             end

                                

           %% Salary Increase
            obj.StartSal = 50000;
            obj.meric = 0.005;
            obj.HistInf = 0.02;
            
            %Set up salary structure at inception, calculate total starting payroll
            obj.salary0 = zeros(obj.w, 1);
            obj.salary0(obj.EntryAge) = obj.StartSal;
            obj.salary0((obj.EntryAge + 1):(obj.RetAge - 1)) = obj.salary0(obj.EntryAge) * (1 + obj.meric) .^ (1:(obj.RetAge - obj.EntryAge - 1));
            obj.salary0(obj.RetAge:obj.w) = 0;
            obj.StartPayroll = sum(obj.nStart .* obj.salary0);
            
%             obj.CumInf = csvread("../../Data/AssetInput/CumInf.csv");
            %salary[x,t,s] is the salary rate applicable during the following year to a member aged x at time t under scenario s
            obj.CumInf = zeros(obj.ProjYear, obj.nscen);
            for i = 1:obj.nscen
                obj.CumInf(:, i) = cumprod(1 + obj.InfAnnual(:, i));
            end           
            
            obj.salary = reshape(kron(obj.CumInf, obj.salary0), [size(obj.salary0,1), size(obj.CumInf,1), size(obj.CumInf,2)]);
            
           %% Benefit Target Function
            %valuation method to be added
            %asset mix
            obj.TargetCalAssump = TargetCalAssump;
%             obj.StockWeight = TargetCalAssump.StockWeight; %0.5;
%             obj.BondWeight = TargetCalAssump.BondWeight; %0.5;
%             obj.ValRateTypeContri = TargetCalAssump.ValRateType; %"BondBased";
%             obj.EquityRPContri = 0.02;
%             obj.LboundContri = 0.01;
%             obj.UboundContri = 0.1;
%             obj.FixValRateContri = 0.05;      
            
            obj.TargetBen = TargetCalAssump.TargetBen;
            obj.ContriRate = TargetCalAssump.ContriRate;
            
            if(TargetCalAssump.SetContri == true)
                obj.BenTargetFunc(Obj.ContriRate);
            else
                obj.ContriRateFunc(obj.TargetBen);
            end
            
                      
            obj.ValRateFunc();
            obj.AssetReturnFunc(obj.PlanDesign.StockWeight, obj.PlanDesign.BondWeight);
            
            %% IC
            obj.InitialFundFunc();        
            
           %% PUC closed Group (function of scenario)
%             obj.NoActionRange= true;
%             obj.TriggerLFR = 0.9;
%             obj.TriggerUFR = 1.1;
%             obj.AdjAction = "edge";
%             obj.AdjLFR = 0.9;
%             obj.AdjUFR = 1.1;
            

            obj.SimulationResult = struct('BenefitPaidOutBOY', zeros(obj.ProjYear+1, obj.nscen), 'ContribBOY', zeros(obj.ProjYear+1, obj.nscen),...
                    'TotalAssetBOY', zeros(obj.ProjYear+1, obj.nscen), 'TotalAssetEOY', zeros(obj.ProjYear, obj.nscen), ...                                       
                    'BenPaid', zeros(obj.ProjYear+1, obj.ProjYear+1, obj.nscen),...
                    'ContribPaid', zeros(obj.ProjYear+1, obj.ProjYear+1, obj.nscen), 'FundGenBOY', zeros(obj.ProjYear+1, obj.ProjYear+1, obj.nscen),...
                    'FundGenEOY', zeros(obj.ProjYear+1, obj.ProjYear+1, obj.nscen));
                
            for i=1:obj.nscen
                obj.DCProjection(i);
            end    
            
%             obj.ProbIncDecFunc();
%             
            obj.ReplaceRatioFunc();
        end
        
        %% Function Helper
        function [AnnFactorsContri, AnnFactorsSalInc] = BTCRFuncHelper(obj)
            if obj.TargetCalAssump.ValRateType == "BondBased"
                CurrValRate =min( max(exp(obj.z0(2)*12)-1, obj.TargetCalAssump.Lbound), obj.TargetCalAssump.Ubound);
            elseif obj.TargetCalAssump.ValRateType == "FixValRate"
                CurrValRate = obj.TargetCalAssump.FixValRate;
            else	
                CurrValRate = min(max(obj.TargetCalAssump.StockWeight*(exp(obj.z0(2)*12)-1+obj.TargetCalAssump.EquityRP) + obj.TargetCalAssump.BondWeight*(exp(obj.z0(2)*12)-1) + 0.0025, obj.TargetCalAssump.Lbound), obj.TargetCalAssump.Ubound);
            end

            RateIndex = find(obj.rates<=CurrValRate, 1, 'last');
            if isempty(RateIndex)
                RateIndex = 0;
            end
            AnnFactorsContri = obj.AnnFactors(obj.EntryAge, RateIndex);
            RateSalInc = (1 + obj.rates(RateIndex)) / ((1 + obj.meric) * (1 + obj.HistInf)) - 1;
            AnnFactorsSalInc = (1 - (1 / (1 + RateSalInc)) ^ (obj.RetAge - obj.EntryAge)) / (1 - (1 / (1 + RateSalInc)));            
        end
        
        %% Benefit Target Function       
        function obj = BenTargetFunc(obj, ContriRate)           
            [AnnFactorsContri, AnnFactorsSalInc] = obj.BTCRFuncHelper();
            obj.TargetBen = AnnFactorsSalInc * ContriRate / (((1 + obj.meric) * (1 + obj.HistInf)) ^ (obj.RetAge - obj.EntryAge - 1) * (obj.RetAge - obj.EntryAge) * AnnFactorsContri);
        end
        
        %% ContriRate Function
         function obj = ContriRateFunc(obj, TargetBen)           
            [AnnFactorsContri, AnnFactorsSalInc] = obj.BTCRFuncHelper();
            obj.ContriRate = TargetBen * (((1 + obj.meric) * (1 + obj.HistInf)) ^ (obj.RetAge - obj.EntryAge - 1) * (obj.RetAge - obj.EntryAge) * AnnFactorsContri) / AnnFactorsSalInc;
         end      
                
        %% Valuation Rate Setting
       %bond based valuation: use the 15 year bond yield as the valuation rate 
       %Annual compounted bond yeild for 15 year zero coupon bond 
       %The valuation bond yield is floored at 0.01, and the maximum will be 0.1
       function obj = ValRateFunc(obj) 
            if obj.PlanDesign.ValRateType == "BondBased"
                obj.ValRateMat = obj.BondYield15;
                obj.InitialValRate = exp(obj.z0(2)*12)-1;
            elseif obj.PlanDesign.ValRateType == "FixValRate"
                obj.ValRateMat = ones(obj.ProjYear+1, obj.nscen) .* obj.PlanDesign.FixValRate;
                obj.InitialValRate = obj.PlanDesign.FixValRate;
            else
                obj.ValRateMat = (obj.BondYield15+obj.PlanDesign.EquityRP) .* obj.PlanDesign.StockWeight + obj.BondYield15 .* obj.PlanDesign.BondWeight + 0.0025;
                obj.InitialValRate = (exp(obj.z0(2)*12)-1 + obj.PlanDesign.EquityRP)* obj.PlanDesign.StockWeight +(exp(obj.z0(2)*12)-1)* obj.PlanDesign.BondWeight + 0.0025;
            end

            obj.ValRateMat(obj.ValRateMat < obj.PlanDesign.Lbound) = obj.PlanDesign.Lbound;
            obj.ValRateMat(obj.ValRateMat > obj.PlanDesign.Ubound) = obj.PlanDesign.Ubound;
            
            obj.InitialValRate = min(max(obj.InitialValRate, obj.PlanDesign.Lbound),obj.PlanDesign.Ubound);
       end
        
        %% AssetReturnFunc
        function obj = AssetReturnFunc(obj, StockWeight, BondWeight)
            %obj.ARMat = (obj.SrAnnual + obj.DivAnnual) .* StockWeight + obj.BondAnnual .* BondWeight;
            obj.ARMat = (obj.SrAnnual) .* StockWeight + obj.BondAnnual .* BondWeight;
        end
        
        %% InitialLia Inidividual Contribution Approach
        function obj = InitialFundFunc(obj)
            Salary0 = obj.salary0(obj.EntryAge:obj.w);
%            RateIndex = find(obj.rates<=obj.InitialValRate, 1, 'last');
            
            if obj.TargetCalAssump.ValRateType == "BondBased"
                CurrValRate =min( max(exp(obj.z0(2)*12)-1, obj.TargetCalAssump.Lbound), obj.TargetCalAssump.Ubound);
            elseif obj.TargetCalAssump.ValRateType == "FixValRate"
                CurrValRate = obj.TargetCalAssump.FixValRate;
            else
                CurrValRate = min(max(obj.TargetCalAssump.StockWeight*(exp(obj.z0(2)*12)-1+obj.TargetCalAssump.EquityRP) + obj.TargetCalAssump.BondWeight*(exp(obj.z0(2)*12)-1) + 0.0025, obj.TargetCalAssump.Lbound), obj.TargetCalAssump.Ubound);
            end
            
            RateIndex = find(obj.rates<=CurrValRate, 1, 'last');
            
            AnnFactor = obj.AnnFactors((obj.EntryAge:obj.w), RateIndex);
            
            ProjFAE1 = zeros(obj.w - obj.EntryAge + 1, 1);
            for i =1:(obj.RetAge - obj.EntryAge)
                ProjFAE1(i) = Salary0(i) * ((1 + obj.meric) * (1 + obj.HistInf)) ^ (obj.RetAge - obj.EntryAge - i);
            end
            
            for i=(obj.RetAge - obj.EntryAge + 1):(obj.w - obj.EntryAge + 1)
                ProjFAE1(i) = Salary0(obj.RetAge - obj.EntryAge) * (1 + obj.HistInf) ^ (-i + (obj.RetAge - obj.EntryAge));
            end   
            
            %PVB + PVFB
            
            PVFB = obj.TargetBen .* ProjFAE1 *35 .* AnnFactor .* obj.nStart(obj.EntryAge:obj.w);
            
            %Annuity Factor for discounting Future Contribution taking into
            %account of the salary increase and the discount rate is
            %consistent with the one discounting benefit
            RateSalInc =  (1 + obj.rates(RateIndex)) / ((1 + obj.meric) * (1 + obj.HistInf)) - 1;
            DisctSal = 1 / (1 + RateSalInc);           
            Index = (obj.RetAge - obj.EntryAge):-1:1;
            AnnFactorContri = zeros(obj.w-obj.EntryAge+1,1);
            AnnFactorContri(1:(obj.RetAge-obj.EntryAge),1) = (1-(DisctSal.^Index))/(1 - DisctSal);
            
            %PVFC 
            PVFC = obj.ContriRate .* Salary0 .* AnnFactorContri .* obj.nStart(obj.EntryAge:obj.w);
            disp(size(PVFC));
            disp(size(PVFB));
            
            obj.InitialFundGen = PVFB - PVFC;
            obj.InitialFund = sum(obj.InitialFundGen);
           
        end
               
        
        
        %% PUCProjection
        function obj = DCProjection(obj, scen)
            SalaryMatrix = obj.salary(obj.EntryAge:obj.w,:,scen);
%            YearOfService = (obj.RetAge - obj.EntryAge) .* ones(obj.w - obj.EntryAge + 1, 1);
%            YearOfService(1:(obj.RetAge - obj.EntryAge)) = 0:(obj.RetAge - obj.EntryAge - 1);
            
            TotalAssetEOY = zeros(obj.ProjYear, 1);
            [ContribBOY, BenefitPaidOutBOY, TotalAssetBOY]= deal(zeros(obj.ProjYear+1, 1));

            [BenPaid, AnnFacMat, FundGenBOY, FundGenEOY] = deal(zeros(obj.w - obj.EntryAge + 1, obj.ProjYear + 1));
            ContribPaid = zeros(obj.w - obj.EntryAge + 1, obj.ProjYear+1);

            AssetReturn = obj.ARMat(:,scen);
            %Set Up Annuity Factor
            ValRate = obj.ValRateMat(:,scen);

            for i = 1:(obj.ProjYear + 1)
                RateIndex = find(obj.rates<=ValRate(i), 1, 'last');
                if isempty(RateIndex)
                    RateIndex = 0;
                end
                AnnFacMat(:, i) = obj.AnnFactors((obj.EntryAge:obj.w), RateIndex);
            end

            %Set Up Contribution Vector
            Salary0 = obj.salary0(obj.EntryAge:obj.w);
            ContribBOY(1)= obj.ContriRate * obj.StartPayroll;
            for i = 2:(obj.ProjYear+1)
                ContribBOY(i) = obj.ContriRate * obj.nStart(obj.EntryAge:obj.w)' * SalaryMatrix(:,i-1);
            end
            FundGenEOY(: , 1) = obj.InitialFundGen ./ obj.nStart(obj.EntryAge:obj.w);
            
            ContribPaid(:, 1) = obj.ContriRate * Salary0;
            ContribPaid(:,2:(obj.ProjYear+1)) = obj.ContriRate * SalaryMatrix(:, 1:(obj.ProjYear));


            BenPaid ((obj.RetAge - obj.EntryAge +1):(obj.w-obj.EntryAge+1), 1) ...
                = FundGenEOY((obj.RetAge - obj.EntryAge +1):(obj.w-obj.EntryAge+1), 1)...
                ./AnnFacMat((obj.RetAge - obj.EntryAge +1):(obj.w-obj.EntryAge+1),1);
            
            BenefitPaidOutBOY(1) = sum(BenPaid ((obj.RetAge - obj.EntryAge +1):(obj.w-obj.EntryAge+1), 1) .* obj.nStart(obj.RetAge:obj.w));

            FundGenBOY(: , 1) = FundGenEOY (: ,1) - BenPaid(: ,1) + ContribPaid(: , 1);
                     
            TotalAssetBOY(1) = sum(FundGenBOY(:,1).* obj.nStart(obj.EntryAge:obj.w));

            for i = 2:(obj.ProjYear + 1)
  
                TotalAssetEOY(i - 1) = TotalAssetBOY(i - 1) .* (1+AssetReturn(i - 1));
                FundGenEOY(1 , i) = 0;
                FundGenEOY(2:(obj.w-obj.EntryAge+1), i) = FundGenBOY(1:(obj.w-obj.EntryAge),i-1) .*(1+AssetReturn(i - 1));  
                
                %cautious about the last year: HAVE LAST YEAR PAYMENT
                BenPaid ((obj.RetAge - obj.EntryAge +1):(obj.w-obj.EntryAge+1), i)...
                    = FundGenEOY((obj.RetAge - obj.EntryAge +1):(obj.w-obj.EntryAge+1), i)...
                    ./AnnFacMat((obj.RetAge - obj.EntryAge +1):(obj.w-obj.EntryAge+1),i);
                
                FundGenBOY(:, i) = FundGenEOY(:,i) - BenPaid(:,i) + ContribPaid(:,i);

                BenefitPaidOutBOY(i) = sum(BenPaid ((obj.RetAge - obj.EntryAge +1):(obj.w-obj.EntryAge+1), 1) .* obj.nStart(obj.RetAge:obj.w));
                TotalAssetBOY(i) = sum(FundGenBOY(:,i).* obj.nStart(obj.EntryAge:obj.w));
            end


            obj.SimulationResult.BenefitPaidOutBOY(:,scen) = BenefitPaidOutBOY;
            obj.SimulationResult.BenPaid(:,:,scen) = BenPaid;
            obj.SimulationResult.ContribBOY(:,scen) = ContribBOY;
            obj.SimulationResult.TotalAssetBOY(:,scen) = TotalAssetBOY;
            obj.SimulationResult.TotalAssetEOY(:,scen) = TotalAssetEOY;
            obj.SimulationResult.ContribPaid(:,:,scen) = ContribPaid;
            obj.SimulationResult.FundGenBOY(:,:,scen) =FundGenBOY;           
            obj.SimulationResult.FundGenEOY(:,:,scen) =FundGenEOY;
        end
        
        %% ProbIncDecFunc
        function obj = ProbIncDecFunc(obj)
            ProbIncDec = zeros(3, obj.ProjYear);
            FRMat = zeros(obj.nscen, obj.ProjYear);
            DistBenAdj = zeros(9, obj.ProjYear);
            %AdjustMat = zeros(obj.nscen, obj.ProjYear);

            for i = 1:obj.nscen
                FRMat(i,:) = obj.SimulationResult.FundRatio(:,i);
            end

            if (obj.PlanDesign.NoActionRange) 
                for i = 1:obj.ProjYear
                    ProbIncDec(1, i) = sum(FRMat(:, i) < obj.PlanDesign.TriggerLFR) / obj.nscen;
                    ProbIncDec(2, i) = sum(FRMat(:, i) >= obj.PlanDesign.TriggerLFR & FRMat(:, i) <= obj.PlanDesign.TriggerUFR) / obj.nscen;
                    ProbIncDec(3, i) = sum(FRMat(:, i) > obj.PlanDesign.TriggerUFR) / obj.nscen;
                end
            else
                for i = 1:obj.ProjYear
                    ProbIncDec(1, i) = sum(FRMat(:, i) < 1) / obj.nscen;
                    ProbIncDec(2, i) = sum(FRMat(:, i) == 1) / obj.nscen;
                    ProbIncDec(3, i) = sum(FRMat(:, i) > 1) / obj.nscen;
                end
            end

            if obj.PlanDesign.NoActionRange
                AdjustMat = FRMat;
                if obj.PlanDesign.AdjAction ~= 'half'
                    AdjustMat(FRMat > obj.PlanDesign.TriggerUFR) = AdjustMat(FRMat > obj.PlanDesign.TriggerUFR) ./ obj.PlanDesign.AdjUFR - 1;
                else
                    obj.PlanDesign.AdjUFR = 0.5 * (obj.PlanDesign.TriggerUFR + FRMat(FRMat > obj.PlanDesign.TriggerUFR));
                    AdjustMat(FRMat > obj.PlanDesign.TriggerUFR) = AdjustMat(FRMat > obj.PlanDesign.TriggerUFR) ./ obj.PlanDesign.AdjUFR - 1;
                end

                AdjustMat(FRMat < obj.PlanDesign.TriggerLFR) = AdjustMat(FRMat < obj.PlanDesign.TriggerLFR) ./ obj.PlanDesign.AdjLFR - 1;
                AdjustMat(FRMat >= obj.PlanDesign.TriggerLFR & FRMat <= obj.PlanDesign.TriggerUFR) = 0;
            else
                AdjustMat = FRMat;
                AdjustMat = AdjustMat-1;
            end

            for i = 1:obj.ProjYear
                DistBenAdj(1, i) = sum(AdjustMat(:, i) < -0.1) / obj.nscen;
                DistBenAdj(2, i) = sum(AdjustMat(:, i) >= -0.1 & AdjustMat(:, i) < -0.05) / obj.nscen;
                DistBenAdj(3, i) = sum(AdjustMat(:, i) >= -0.05 & AdjustMat(:, i) < -0.02) / obj.nscen;
                DistBenAdj(4, i) = sum(AdjustMat(:, i) >= -0.02 & AdjustMat(:, i) < 0) / obj.nscen;
                DistBenAdj(5, i) = sum(AdjustMat(:, i) == 0) / obj.nscen;
                DistBenAdj(6, i) = sum(AdjustMat(:, i) > 0 & AdjustMat(:, i) < 0.02) / obj.nscen;
                DistBenAdj(7, i) = sum(AdjustMat(:, i) >= 0.02 & AdjustMat(:, i) < 0.05) / obj.nscen;
                DistBenAdj(8, i) = sum(AdjustMat(:, i) >= 0.05 & AdjustMat(:, i) < 0.1) / obj.nscen;
                DistBenAdj(9, i) = sum(AdjustMat(:, i) >= 0.1) / obj.nscen;
            end

            obj.DistReturnList = struct('ProbIncDec',ProbIncDec,'DistBenAdj',DistBenAdj);            
        end
        
        %% ReplaceRatioFunc
        function obj = ReplaceRatioFunc(obj)
            PensionAtRetire = zeros(obj.nscen, obj.ProjYear);
            FinalSalary = zeros(obj.nscen, obj.ProjYear);
            %ReplaceRatioAtRetire = zeros(obj.nscen, obj.ProjYear);
            for i=1:obj.nscen
                PensionAtRetire(i,:) = obj.SimulationResult.BenPaid(36, 2:(obj.ProjYear+1), i);
                FinalSalary(i,:) = obj.salary(64,:,i);
            end
            obj.ReplaceRatioAtRetire = PensionAtRetire ./ FinalSalary;
        end
    end
end
