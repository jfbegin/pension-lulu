
classdef SimResultRN < handle
%{
Class for RN estimation

function yStruct = termstrt(obj, t, lambda)
    Normally outputs the bond yield of different maturities (i.e. 12 36 60
    84 120 144 168 180 months) at time t
    However, for testing reasons, it will output Z_t when NaN value occurs
    and interrupt the function.
    
function outMat = test(obj, t, lambda, ts, scn) 
    Outputs Z_t and Sigma_t for the scenario scn up to time step ts
    For testing purpose only
    
function yTermLSE = RNEst(obj, lambda)
    Outputs the final squared difference between the model results and the
    historical data
    takes around 60 seconds to run

function yTermLSE = RNEstGPU(obj, lambda)
    GPU version of the above function
    takes around 10 seconds to run 
%}
    properties      
        mu        %4*1 vector, the constant term in the VAR equation
        Beta      %4*4 matrix, the auto correlation coefficient
        Cmat      %4*4 diagonal matrix, the constant term in the GARCH equation
        A1        %4*4 diagonal matrix, the ARCH coefficient in the GARCH equation
        B1        %4*4 diagonal matrix, the GARCH coefficient in the GARCH equation
        LambdaT   %[0, 0, 0, lambda4 - 0.5]
        Sigma0    %Longterm variance
        revlevel  %mean reverting level of the VAR model
        T         %maturities in month
        L4Vec     
        LambdaRN  
        lambda0    % -- new
        lambda1    % -- new
        
        simZt      % -- new
        simHt      % -- new
        simRest    % -- new
        
        termStruct % -- new
        numScen    % -- new
        simIdx     % -- new
        
        % inovs and HistData are passed by reference
        inovs     %inovation terms generate from Normal(0, 1), 
                  %4 dimensions:10000 scenarios for 4 variables for 180 months starting from 301 time points in historical data
        HistData  %Historical Data
        HistRest  %Historical Residual 
        HistHMat  %Historical Ht Matrix
        Sigma     %Sigma from VAR estimation
        HistTStrtData %Historical Data for bond yield at different maturities
        
        boundedMat
        boundedMatScen
    end

    
    methods        
        function obj = SimResultRN(suminfo, HistData, Sigma, revlevel, inovs, HistTStrtData, lambda)
            % Constructor
            if nargin > 0
                obj.numScen = 10000;
                obj.revlevel = revlevel;
                obj.update(suminfo, lambda);
                obj.HistData = HistData;
                obj.Sigma = Sigma;
                obj.inovs = inovs;
                obj.T = [12 24 36 48 60 72 84 96 108 120 132 144 156 168 180];
                obj.HistHMat = suminfo.CondVariances;
                obj.HistRest = suminfo.Residuals;
                obj.L4Vec =  [0; 0; 0; suminfo.lambda4];
                %obj.LambdaRN = diag(obj.LambdaT - obj.L4Vec);
                obj.LambdaRN = [0; 0; 0; -0.5];
                obj.HistTStrtData = HistTStrtData(:, obj.T / 12);
                
                %obj.boundedMat = gpuArray(25.*repmat(obj.Sigma0,661,obj.numScen));
                obj.boundedMat = gpuArray(25.*repmat(obj.Sigma0,56,obj.numScen));
                obj.boundedMatScen = 25.*repmat(obj.Sigma0, 1, 10000);
                
%                obj.simZt = zeros(662, 5, 10000);
%                obj.simRest = zeros(661, 5, 10000);
%                obj.simHt = zeros(661, 5, 10000);
                
                obj.simZt = zeros(663, 4, 10000);
                obj.simRest = zeros(662, 4, 10000);
                obj.simHt = zeros(662, 4, 10000);
                
                %obj.termStruct = zeros(661, 15, obj.numScen);
                obj.termStruct = zeros(56, 15, obj.numScen);
                obj.simIdx = 0;
            end
        end

    
        function [mu, Beta, Cmat, A1, B1, LambdaT, Sigma0, lambda0, lambda1] = load(obj, suminfo, lambda)
            % Load parameter vector into corresponding variables
            Beta = reshape(suminfo.xCenter(1:16),[4,4]);
            Cmat = diag(suminfo.xCenter(17:20));
            A1 = diag(suminfo.xCenter(21:24));
            B1 = diag(suminfo.xCenter(25:28));
            
            LambdaT = [0; 0; 0; suminfo.lambda4-0.5];
            
            mu = (eye(4) - Beta) * obj.revlevel;
            
            Sigma0 = suminfo.xCenter(17:20)./(ones(4,1)-suminfo.xCenter(21:24)-suminfo.xCenter(25:28));
            
            lambda0 = [lambda(1); lambda(2); 0; mu(4)]; % modified
            lambda1 = [lambda(3:6);lambda(7:10);0 0 0 0; Beta(4,:)]; % modified
        end

        function obj = update(obj, suminfo, lambda)
            % Update parameters
            [obj.mu, obj.Beta, obj.Cmat, obj.A1, obj.B1, obj.LambdaT, obj.Sigma0, obj.lambda0, obj.lambda1] = obj.load(suminfo, lambda);
        end
        
        function obj = genScenario(obj)
%            Z0 = obj.HistData.m(301,:)';
%            startz0 = [-8.0293;-6.1287;0.0024; 0.003; -7.0583];
            startz0 = obj.revlevel;
            Z0 = startz0;
            
            obj.simZt(1,:,:) = repmat(Z0, 1, 10000);
            
%            Z1 = obj.HistData.m(302,:)';
            Z1 = startz0;
            
            rng('default');
%           inovt = randn(5,10000,660,'single');
            inovt = randn(4,10000,661,'single');
            
%            LastResT = obj.HistRest(301,:)';
%            LastSigT = obj.HistHMat(301,:)';
            
            LastResT = [0;0;0;0];
            LastSigT = obj.Sigma0;
                                 
            SigmaT = diag(obj.Cmat) + diag(obj.A1) .* ...
                     ((LastResT - diag(obj.L4Vec) * LastSigT - (obj.lambda0 + obj.lambda1 * Z0)).^2)...
                     + diag(obj.B1) .* LastSigT; 
                 
            SigmaT = repmat(SigmaT, 1, 10000);
            LastSigT = repmat(LastSigT, 1, 10000);
            obj.simHt(1,:,:) = LastSigT;
            %obj.simHt(2,:,:) = SigmaT;
            LastZT = repmat(Z1, 1, 10000);     
            obj.simZt(2,:,:) = LastZT;
            obj.simRest(1,:,:) = repmat(LastResT, 1, 10000);          

            %for i=2:661
            for i=2:662
                obj.simHt(i,:,:) = SigmaT;
                Yt = sqrt(SigmaT) .* inovt(:,:,i-1);   
                obj.simRest(i,:,:) = Yt;
                
                %disp(size(diag(obj.Beta) .* obj.simZt(i-1,:,:)));
                %disp(size(Yt));
                %disp(size(obj.LambdaT .* SigmaT));
             
                Zt =  ((SigmaT .* obj.LambdaRN) + (obj.mu - obj.lambda0)) + ((obj.Beta - obj.lambda1)* squeeze(obj.simZt(i,:,:)) + Yt); % modified
                
                SigmaT = obj.boundedSigmaScen( diag(obj.Cmat) + diag(obj.A1) .* ...
                         ((Yt - SigmaT.*obj.L4Vec - (obj.lambda1 * squeeze(obj.simZt(i,:,:)) +obj.lambda0)).^2) ...
                         + diag(obj.B1) .* SigmaT );
                
                obj.simZt(i+1,:,:) = Zt;
            end
        end        

        function yStruct = oneSceBondYield(obj, scen)
            m_idx = (1:4:2644)';
            A = gpuArray(repmat(diag(obj.A1),661,1));
            B = gpuArray(repmat(diag(obj.B1),661,1));
            C = gpuArray(repmat(diag(obj.Cmat),661,1));
            L4V = gpuArray(repmat(obj.L4Vec,661,1));
            LRN = gpuArray(repmat(obj.LambdaRN, 661, 1));
                   
            Z0 = gpuArray(reshape(squeeze(obj.simZt(1:661,:,scen))', 2644, 1));
            Z1 = gpuArray(reshape(squeeze(obj.simZt(2:662,:,scen))', 2644, 1));
          
            L0 = gpuArray(repmat(obj.lambda0, 661, 1));
            L1 = gpuArray(kron(eye(661), obj.lambda1));
            mu_L0 = gpuArray(repmat(obj.mu - obj.lambda0, 661, 1));
            B_L1 = gpuArray(kron(eye(661), obj.Beta - obj.lambda1));
            
            LastResT = gpuArray(reshape(squeeze(obj.simRest(:,:,scen))', 2644, 1));
            LastSigT = gpuArray(reshape(squeeze(obj.simHt(:,:,scen))', 2644, 1));
            
            SigmaT = gpuArray(C + A .* ((LastResT - L4V .* LastSigT - (L0 + L1 * Z0)).^2)...
                     + B .* LastSigT);
                 
            SigmaT = gpuArray(repmat(SigmaT, 1, obj.numScen));
            LastZT = gpuArray(repmat(Z1, 1, obj.numScen));          
            yStruct = gpuArray(zeros(661, size(obj.T, 2)));
            sumZT = gpuArray(exp(LastZT(m_idx,:)));
            idx = 1;
            
            for i=2:180
                Yt = sqrt(SigmaT) .* obj.inovs.m(:,:,i-1);
                Zt =  (SigmaT .* LRN + mu_L0) + B_L1 * LastZT + Yt;           
                sumZT = sumZT + exp(Zt(m_idx,:));
                if i == obj.T(idx)
                    yStruct(:,idx) = -log(mean(exp(-sumZT),2))/obj.T(idx)*12;
                    idx = idx + 1;
                end

                SigmaT = obj.boundedSigma(C + A .* ...
                         ((Yt - SigmaT .* L4V - (L1 * LastZT + L0)).^2) ...
                         + B .* SigmaT);               
                LastZT = Zt;
            end
            yStruct = gather(yStruct);
            %yTermLSE = gather(sum(sum((obj.HistTStrtData - yStruct).^2)));            
        end

        function yStruct = oneSceBondYieldYear(obj, scen)
            m_idx = (1:4:224)';
            A = gpuArray(repmat(diag(obj.A1),56,1));
            B = gpuArray(repmat(diag(obj.B1),56,1));
            C = gpuArray(repmat(diag(obj.Cmat),56,1));
            L4V = gpuArray(repmat(obj.L4Vec,56,1));
            LRN = gpuArray(repmat(obj.LambdaRN, 56, 1));
                   
%            Z0 = gpuArray(reshape(squeeze(obj.simZt(1:12:661,:,scen))', 280, 1));
%            Z1 = gpuArray(reshape(squeeze(obj.simZt(2:12:662,:,scen))', 280, 1));
            
            Z0 = gpuArray(reshape(squeeze(obj.simZt(2:12:662,:,scen))', 224, 1));
            Z1 = gpuArray(reshape(squeeze(obj.simZt(3:12:663,:,scen))', 224, 1));
          
            L0 = gpuArray(repmat(obj.lambda0, 56, 1));
            L1 = gpuArray(kron(eye(56), obj.lambda1));
            mu_L0 = gpuArray(repmat(obj.mu - obj.lambda0, 56, 1));
            B_L1 = gpuArray(kron(eye(56), obj.Beta - obj.lambda1));
            
%            LastResT = gpuArray(reshape(squeeze(obj.simRest(1:12:661,:,scen))', 280, 1));
%            LastSigT = gpuArray(reshape(squeeze(obj.simHt(1:12:661,:,scen))', 280, 1));
            
            LastResT = gpuArray(reshape(squeeze(obj.simRest(2:12:662,:,scen))', 224, 1));
            LastSigT = gpuArray(reshape(squeeze(obj.simHt(2:12:662,:,scen))', 224, 1));
            
            SigmaT = gpuArray(C + A .* ((LastResT - L4V .* LastSigT - (L0 + L1 * Z0)).^2)...
                     + B .* LastSigT);
                 
            SigmaT = gpuArray(repmat(SigmaT, 1, obj.numScen));
            LastZT = gpuArray(repmat(Z1, 1, obj.numScen));          
            yStruct = gpuArray(zeros(56, size(obj.T, 2)));
            sumZT = gpuArray(LastZT(m_idx,:));
            idx = 1;
            
            for i=2:180
                Yt = sqrt(SigmaT) .* obj.inovs.m(:,:,i-1);
                Zt =  (SigmaT .* LRN + mu_L0) + B_L1 * LastZT + Yt;           
                sumZT = sumZT + Zt(m_idx,:);
                if i == obj.T(idx)
                    yStruct(:,idx) = -log(mean(exp(-sumZT),2))/obj.T(idx)*12;
                    idx = idx + 1;
                end

                SigmaT = obj.boundedSigma(C + A .* ...
                         ((Yt - SigmaT .* L4V - (L1 * LastZT + L0)).^2) ...
                         + B .* SigmaT);               
                LastZT = Zt;
            end
            yStruct = gather(yStruct);
            %yTermLSE = gather(sum(sum((obj.HistTStrtData - yStruct).^2)));            
        end        
        
        function boundedSigmas = boundedSigmaScen(obj, Sigmas)
                %[~,N] = size(Sigmas);
                boundedSigmas = max(eps,min(Sigmas, obj.boundedMatScen));
        end        
        
        function boundedSigmas = boundedSigma(obj, Sigmas)
                %[~,N] = size(Sigmas);
                boundedSigmas = max(eps,min(Sigmas, obj.boundedMat));
        end
        
        function obj = fillTermStruct(obj, numIter)
            if obj.simIdx >= 10000
                return
            end
            
            if obj.simIdx + numIter > 10000 
                numIter = 10000 - obj.simIdx;
                fprintf('Number of iteration: %d\n', numIter);
            end
           
            for i = (obj.simIdx+1):(obj.simIdx+numIter)
                fprintf('iteration: %d\n', i);             
                %obj.termStruct(:,:,i) = obj.oneSceBondYield(i);
                obj.termStruct(:,:,i) = obj.oneSceBondYieldYear(i);
                obj.simIdx = obj.simIdx + 1;
            end
        end
    end
end


