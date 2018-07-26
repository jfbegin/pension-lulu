classdef AssetScenario < handle
    % This class summarizes and aggregates the output from the ESG, making
    % it easier in the liabilities calculations. 
    properties
        Zt
        nScen
        projYear
        termStruct
        
        srAnnualConti  %55*10000 matrix
        %divAnnualConti
        bondAnnualConti
        infAnnualConti
        
        srAnnualEfft %
        %divAnnualEfft %
        bondAnnualEfft%
        infAnnualEfft%
        
        bondYield_S_Efft %56*10000 matrix
        bondYield_L_Efft %
        
        mt
        Mt%
        
    end
    
    methods
        function obj = AssetScenario(Zt, nScen, projYear, termStruct)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Zt = Zt;
            obj.nScen = nScen;
            obj.projYear = projYear;
            obj.termStruct = termStruct;
            
            obj.StockReturn();
            %obj.DivYield();
            obj.BondReturn();
            obj.Inflation();
            obj.BondYield();
            obj.DiscountFactor();
        end
        
        function obj = StockReturn(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            srMonth = squeeze(obj.Zt(:,1,:) + obj.Zt(:,4,:));
%            obj.srAnnualConti = squeeze(sum(reshape(srMonth(3:end,:),12,55,obj.nScen)));
            obj.srAnnualConti = squeeze(sum(reshape(srMonth(3:662,:),12,55,obj.nScen)));
            
            obj.srAnnualEfft = exp(obj.srAnnualConti) - 1;
        end
        
%         function obj = DivYield(obj)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             divMonth = squeeze(obj.Zt(:,5,:));
% %            obj.divAnnualConti = squeeze(sum(reshape(divMonth(3:end,:),12,55,obj.nScen)));
%             obj.divAnnualConti = squeeze(sum(reshape(divMonth(3:662,:),12,55,obj.nScen)));
%             
%             obj.divAnnualEfft = exp(obj.divAnnualConti) - 1;
%         end
        
        function obj = BondReturn(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            bp15t0 = squeeze(obj.termStruct(1:55,15,:));
            bp14t1 = squeeze(obj.termStruct(2:56,14,:));
%             bp15t0 = squeeze(obj.termStruct(1:12:660,15,:));
%             bp14t1 = squeeze(obj.termStruct(13:12:661,14,:));
            obj.bondAnnualConti = bp15t0 * 15 - bp14t1 * 14;            
            obj.bondAnnualEfft = exp(obj.bondAnnualConti) - 1;
        end 
        
        function obj = Inflation(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            infMonth = squeeze(obj.Zt(:,3,:));
%            obj.infAnnualConti = squeeze(sum(reshape(infMonth(3:end,:),12,55,obj.nScen)));
            obj.infAnnualConti = squeeze(sum(reshape(infMonth(3:662,:),12,55,obj.nScen)));
            obj.infAnnualEfft = exp(obj.infAnnualConti) - 1;
        end  
        
        function obj = BondYield(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
%            ZtAnnualStep = obj.Zt(2:12:end,:,:);            
            ZtAnnualStep = obj.Zt(3:12:end,:,:);
            
            obj.bondYield_S_Efft = exp(squeeze(ZtAnnualStep(:,1,:)).*12) - 1;
            obj.bondYield_L_Efft = exp(squeeze(ZtAnnualStep(:,2,:)).*12) - 1;
        end  
        
        function obj = DiscountFactor(obj)
            obj.mt = zeros(56, obj.nScen);
            temp = squeeze(obj.Zt(:,1,:));
%            obj.mt(2:end, :) = squeeze(sum(reshape(temp(2:661,:),12,55,obj.nScen)));
            obj.mt(2:end, :) = squeeze(sum(reshape(temp(3:662,:),12,55,obj.nScen)));
            
            obj.Mt = cumsum(obj.mt);
            obj.Mt = exp(-obj.Mt);
        end
    end
end

