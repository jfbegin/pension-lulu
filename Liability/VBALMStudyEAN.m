classdef VBALMStudyEAN < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ALMResult
        DCALMResult
        
        CashflowBasic
        ResiValueAtEnd
        InitialContriPaid
        GenAccountScen
        VGenAccountBasic
        VGenAccountRes
        VGenAccountBen
        
        DCCashflowBasic
        DCResiValueAtEnd
        DCInitialContriPaid
        
        DCGenAccountScen
        DCGenResidualScen
        DCGenBenScen
        
        VDCGenAccountRes
        VDCGenAccountBen
        VDCGenAccountBasic
        
        SurplusOptMat
        DeficitOptMat
        ResSurplus
        ResDeficit
        
        GenSurplusScen
        GenDeficitScen
        GenResSurplusScen
        GenResDeficitScen
        GenTotalSurplusScen
        GenTotalDeficitScen
        
        GenSurplus
        GenDeficit
        GenResSurplus
        GenResDeficit
        GenTotalSurplus
        GenTotalDeficit
        
        TermEndAccLia       
        
        InitialAccLia
        
        GenResidualScen
        GenBenScen
    end
    
    methods
        function obj = VBALMStudyEAN(ALMResult ,DCALMResult)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.ALMResult = ALMResult;
            obj.DCALMResult = DCALMResult;
            obj.GenerationCashFlow();
            obj.DCGenerationCashFlow();
            obj.GenerationValueFunc();
            obj.DCGenerationValueFunc();
            obj.SurplusDeficitFunc();            
        end
        
        function obj = GenerationCashFlow(obj)
            w = obj.ALMResult.w;
            EntryAge = obj.ALMResult.EntryAge;
            ProjYear = obj.ALMResult.ProjYear;
            RetAge = obj.ALMResult.RetAge;
            nscen = obj.ALMResult.nscen;
            
            disp("run Gen Func");
%            obj.CashflowBasic = zeros (w + (ProjYear - EntryAge), ProjYear, nscen);
            obj.CashflowBasic = zeros (w + (ProjYear - EntryAge) +1 , ProjYear+1, nscen);
            for i=1:nscen
%                disp(i);
%                CashflowMat = zeros(w-EntryAge+1, ProjYear);
                CashflowMat = zeros(w-EntryAge+1, ProjYear+1);

%                CashflowMat(1:(RetAge - EntryAge),:) = -1 * ...
%                   obj.ALMResult.SimulationResult.ContribPaid(1:(RetAge -
%                   EntryAge),(1:ProjYear),i);
                CashflowMat(1:(RetAge - EntryAge),:) = -1 * ...
                    obj.ALMResult.SimulationResult.ContribPaid(1:(RetAge - EntryAge),:,i);
                
%                CashflowMat((RetAge - EntryAge + 1):(w - EntryAge + 1),:) = ...
%                    obj.ALMResult.SimulationResult.AdjAccBen((RetAge - EntryAge + 1):(w - EntryAge + 1), (1:ProjYear),i)
                CashflowMat((RetAge - EntryAge + 1):(w - EntryAge + 1),:) = ...
                    obj.ALMResult.SimulationResult.AdjAccBen((RetAge - EntryAge + 1):(w - EntryAge + 1), :,i);
                
%                for j=1:ProjYear
                for j=1:(ProjYear+1)
%                    obj.CashflowBasic((ProjYear - j + 1):(w + ProjYear - EntryAge - j + 1), j, i) = ...
%                        CashflowMat(:, j);
                    obj.CashflowBasic((ProjYear+1 - j + 1):(w + ProjYear+1 - EntryAge - j + 1), j, i) = ...
                        CashflowMat(:, j);
                end             
            end
            
            %% CHANGE WITH PUC !
            % The only difference with VBALMStudy: in EAN, we have some 
            % estimator of the residual amount, which we did not have with 
            % the PUC method.
            
            obj.TermEndAccLia = zeros(w + (ProjYear - EntryAge) +1, nscen);
            obj.TermEndAccLia(1:(obj.ALMResult.w-EntryAge +1),(1:nscen)) = obj.ALMResult.SimulationResult.EndLiaInd(:,(1:nscen));
            % the contribution is negative cash flow, so use negative value
            % to calculate the cum Fund

            for i=1:nscen
%               obj.ResiValueAtEnd(:, i) = obj.ALMResult.SimulationResult.TotalAssetEOY(ProjYear, i) ./ ...
%                   sum(obj.TermEndAccLia(:, i)) .* obj.TermEndAccLia(:, i) ./ 100;
               obj.ResiValueAtEnd(:, i) = obj.ALMResult.SimulationResult.TotalAssetBOY(ProjYear+1, i) ./ ...
                   sum(obj.TermEndAccLia(:, i)) .* obj.TermEndAccLia(:, i) ./ 100;
            end            
            
            %Line Up with DC InitialContri 
            %We are assuming that we are transit from a DC plan to TBP
            % if not doing this the Initial Fund should also be separated
            % based on the DC liability proportionally
%            obj.InitialContriPaid = zeros(w + (ProjYear - EntryAge), nscen);
            obj.InitialContriPaid = zeros(w + (ProjYear - EntryAge)+1, nscen);
            
            for i=1:nscen
               obj.InitialContriPaid((ProjYear+1):(w + ProjYear - EntryAge+1), i) =...
                   obj.ALMResult.InitialLiaICGen ./ obj.ALMResult.nStart(obj.ALMResult.EntryAge:obj.ALMResult.w);             
            end            
        end
        
        % EXACTLY THE SAME AS FOR PUC !
        function obj = DCGenerationCashFlow(obj)
            w = obj.DCALMResult.w;
            EntryAge = obj.DCALMResult.EntryAge;
            ProjYear = obj.DCALMResult.ProjYear;
            RetAge = obj.DCALMResult.RetAge;
            nscen = obj.DCALMResult.nscen;
            
%            obj.DCCashflowBasic = zeros( w + (ProjYear - EntryAge), ProjYear, nscen);
            obj.DCCashflowBasic = zeros( w + (ProjYear - EntryAge)+1, ProjYear, nscen);
            for i=1:nscen
%                DCCashflowMat = zeros(w-EntryAge+1, ProjYear);
                DCCashflowMat = zeros(w-EntryAge+1, ProjYear+1);

%                DCCashflowMat(1:(RetAge - EntryAge),:) = -1 * ...
%                    obj.DCALMResult.SimulationResult.ContribPaid(1:(RetAge - EntryAge),1:ProjYear,i);
                DCCashflowMat(1:(RetAge - EntryAge),:) = -1 * ...
                    obj.DCALMResult.SimulationResult.ContribPaid(1:(RetAge - EntryAge),:,i);
                
%                DCCashflowMat((RetAge - EntryAge + 1):(w - EntryAge + 1),:) = ...
 %                   obj.DCALMResult.SimulationResult.BenPaid((RetAge - EntryAge + 1):(w - EntryAge + 1), (1:ProjYear),i);
                DCCashflowMat((RetAge - EntryAge + 1):(w - EntryAge + 1),:) = ...
                     obj.DCALMResult.SimulationResult.BenPaid((RetAge - EntryAge + 1):(w - EntryAge + 1), :,i);
 
%                    for j=1:(ProjYear) 
                for j=1:(ProjYear+1)
%                    obj.DCCashflowBasic((ProjYear - j + 1):(w + ProjYear - EntryAge - j + 1), j, i) = ...
%                        DCCashflowMat(:, j);
                    obj.DCCashflowBasic((ProjYear+1 - j + 1):(w + ProjYear+1 - EntryAge - j + 1), j, i) = ...
                        DCCashflowMat(:, j);
                end             
            end

%            obj.DCResiValueAtEnd = zeros(w + (ProjYear - EntryAge), nscen);
            obj.DCResiValueAtEnd = zeros(w + (ProjYear - EntryAge)+1, nscen);
            
            for i=1:nscen
%               obj.DCResiValueAtEnd(1:(w - EntryAge), i) = ...
%                   obj.DCALMResult.SimulationResult.FundGenEOY(2:(w - EntryAge + 1), ProjYear + 1, i);
               obj.DCResiValueAtEnd(1:(w - EntryAge), i) = ...
                   obj.DCALMResult.SimulationResult.FundGenBOY(1:(w - EntryAge), ProjYear + 1, i);
            end            
            
%            obj.DCInitialContriPaid = zeros(w + (ProjYear - EntryAge), nscen);
            obj.DCInitialContriPaid = zeros(w + (ProjYear - EntryAge)+1, nscen);
            
            for i=1:nscen
               obj.DCInitialContriPaid((ProjYear+1):(w + ProjYear - EntryAge+1), i) =...
                   obj.DCALMResult.SimulationResult.FundGenEOY(:, 1, i);    
            end            
        end       
        
        
        function obj = GenerationValueFunc(obj)
            w = obj.ALMResult.w;
            EntryAge = obj.ALMResult.EntryAge;
            ProjYear = obj.ALMResult.ProjYear;
            RetAge = obj.ALMResult.RetAge;
            nscen = obj.ALMResult.nscen;          

            
            CumM = obj.ALMResult.Mt;
%            obj.GenAccountScen = zeros(w + (ProjYear - EntryAge), nscen);
%            obj.GenResidualScen = zeros(w + (ProjYear - EntryAge), nscen);
%            obj.GenBenScen = zeros(w+(ProjYear-EntryAge),nscen);
            
            obj.GenAccountScen = zeros(w + (ProjYear - EntryAge)+1, nscen);
            obj.GenResidualScen = zeros(w + (ProjYear - EntryAge)+1, nscen);
            obj.GenBenScen = zeros(w+(ProjYear-EntryAge)+1,nscen);
            
            for i=1:nscen
%               for j=1:(w + (ProjYear - EntryAge))
               for j=1:(w + (ProjYear - EntryAge)+1)
                   %obj.GenResidualScen(j, i) = obj.ResiValueAtEnd(j, i) * CumM(ProjYear + 1, i);
                   obj.GenResidualScen(j, i) = obj.ResiValueAtEnd(j, i) * CumM(ProjYear + 1, i);
                   %obj.GenBenScen(j, i) = sum(obj.CashflowBasic(j,:,i)' .* CumM(1:ProjYear, i));
                   obj.GenBenScen(j, i) = sum(obj.CashflowBasic(j,:,i)' .* CumM(1:(ProjYear+1), i));
                   obj.GenAccountScen(j, i) = obj.GenBenScen(j,i) + obj.GenResidualScen(j,i);
               end
            end
            
            obj.VGenAccountRes = mean(obj.GenResidualScen, 2);
            obj.VGenAccountBen = mean(obj.GenBenScen, 2);
            obj.VGenAccountBasic = mean(obj.GenAccountScen, 2);           
        end
        
        function obj = DCGenerationValueFunc(obj)
            w = obj.ALMResult.w;
            EntryAge = obj.ALMResult.EntryAge;
            ProjYear = obj.ALMResult.ProjYear;
            RetAge = obj.ALMResult.RetAge;
            nscen = obj.ALMResult.nscen;          

            
            CumM = obj.ALMResult.Mt;
%            obj.DCGenAccountScen = zeros(w + (ProjYear - EntryAge), nscen);
%            obj.DCGenResidualScen = zeros(w + (ProjYear - EntryAge), nscen);
%            obj.DCGenBenScen = zeros(w+(ProjYear-EntryAge),nscen);
            
            obj.DCGenAccountScen = zeros(w + (ProjYear - EntryAge)+1, nscen);
            obj.DCGenResidualScen = zeros(w + (ProjYear - EntryAge)+1, nscen);
            obj.DCGenBenScen = zeros(w+(ProjYear-EntryAge)+1,nscen);
            
            for i=1:nscen
%               for j=1:(w + (ProjYear - EntryAge))
               for j=1:(w + (ProjYear - EntryAge)+1)
                   %obj.GenResidualScen(j, i) = obj.ResiValueAtEnd(j, i) * CumM(ProjYear + 1, i);
                   obj.DCGenResidualScen(j, i) = obj.DCResiValueAtEnd(j, i) * CumM(ProjYear + 1, i);
                   %obj.GenBenScen(j, i) = sum(obj.CashflowBasic(j,:,i)' .* CumM(1:ProjYear, i));
                   obj.DCGenBenScen(j, i) = sum(obj.DCCashflowBasic(j,:,i)' .* CumM(1:(ProjYear+1), i));
                   obj.DCGenAccountScen(j, i) = obj.DCGenBenScen(j,i) + obj.DCGenResidualScen(j,i);
               end
            end
            
            obj.VDCGenAccountRes = mean(obj.DCGenResidualScen, 2);
            obj.VDCGenAccountBen = mean(obj.DCGenBenScen, 2);
            obj.VDCGenAccountBasic = mean(obj.DCGenAccountScen, 2);           
        end
        
        function obj = SurplusDeficitFunc(obj)            
            w = obj.ALMResult.w;
            EntryAge = obj.ALMResult.EntryAge;
            ProjYear = obj.ALMResult.ProjYear;
            RetAge = obj.ALMResult.RetAge;
            nscen = obj.ALMResult.nscen; 
            
  %          obj.SurplusOptMat = max((obj.CashflowBasic - obj.DCCashflowBasic), zeros(w + (ProjYear - EntryAge), ProjYear, nscen));
  %          obj.DeficitOptMat = max((obj.DCCashflowBasic - obj.CashflowBasic), zeros(w + (ProjYear - EntryAge), ProjYear, nscen));
  %          obj.ResSurplus = max((obj.ResiValueAtEnd - obj.DCResiValueAtEnd), zeros(w + (ProjYear - EntryAge), nscen));
  %          obj.ResDeficit = max((obj.DCResiValueAtEnd - obj.ResiValueAtEnd), zeros(w + (ProjYear - EntryAge), nscen));
            
            obj.SurplusOptMat = max((obj.CashflowBasic - obj.DCCashflowBasic), zeros(w + (ProjYear - EntryAge)+1, ProjYear+1, nscen));
            obj.DeficitOptMat = max((obj.DCCashflowBasic - obj.CashflowBasic), zeros(w + (ProjYear - EntryAge)+1, ProjYear+1, nscen));
            obj.ResSurplus = max((obj.ResiValueAtEnd - obj.DCResiValueAtEnd), zeros(w + (ProjYear - EntryAge)+1, nscen));
            obj.ResDeficit = max((obj.DCResiValueAtEnd - obj.ResiValueAtEnd), zeros(w + (ProjYear - EntryAge)+1, nscen));
            
            CumM = obj.ALMResult.Mt;
%            obj.GenSurplusScen = zeros(w + (ProjYear - EntryAge), nscen);
%            obj.GenDeficitScen = zeros(w + (ProjYear - EntryAge), nscen);
%            obj.GenResSurplusScen = zeros(w+(ProjYear-EntryAge),nscen);
%            obj.GenResDeficitScen = zeros(w+(ProjYear-EntryAge),nscen);
%            obj.GenTotalSurplusScen = zeros(w+(ProjYear-EntryAge),nscen);
%            obj.GenTotalDeficitScen = zeros(w+(ProjYear-EntryAge),nscen);
            
            obj.GenSurplusScen = zeros(w + (ProjYear - EntryAge)+1, nscen);
            obj.GenDeficitScen = zeros(w + (ProjYear - EntryAge)+1, nscen);
            obj.GenResSurplusScen = zeros(w+(ProjYear-EntryAge)+1,nscen);
            obj.GenResDeficitScen = zeros(w+(ProjYear-EntryAge)+1,nscen);
            obj.GenTotalSurplusScen = zeros(w+(ProjYear-EntryAge)+1,nscen);
            obj.GenTotalDeficitScen = zeros(w+(ProjYear-EntryAge)+1,nscen);
            
            for i=1:nscen
  %             for j=1:(w + (ProjYear - EntryAge))
                for j=1:(w + (ProjYear - EntryAge)+1)
                   obj.GenResSurplusScen(j, i) = obj.ResSurplus(j, i) * CumM(ProjYear + 1, i);
                   obj.GenResDeficitScen(j, i) = obj.ResDeficit(j, i) * CumM(ProjYear + 1, i);
                   obj.GenSurplusScen(j, i) = sum(obj.SurplusOptMat(j,:,i)' .* CumM(1:(ProjYear+1), i));
                   obj.GenDeficitScen(j, i) = sum(obj.DeficitOptMat(j,:,i)' .* CumM(1:(ProjYear+1), i));
                   obj.GenTotalSurplusScen(j,i) = obj.GenSurplusScen(j,i) + obj.GenResSurplusScen(j,i);
                   obj.GenTotalDeficitScen(j,i) = obj.GenDeficitScen(j,i) + obj.GenResDeficitScen(j,i);
               end
            end
            
            
            obj.GenSurplus = mean(obj.GenSurplusScen, 2);
            obj.GenDeficit = mean(obj.GenDeficitScen, 2);
            obj.GenResSurplus = mean(obj.GenResSurplusScen, 2);
            obj.GenResDeficit = mean(obj.GenResDeficitScen, 2);
            obj.GenTotalSurplus = mean(obj.GenTotalSurplusScen, 2);
            obj.GenTotalDeficit = mean(obj.GenTotalDeficitScen, 2);
        end        
        
    end
end

