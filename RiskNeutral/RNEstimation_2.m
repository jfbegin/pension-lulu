
classdef RNEstimation_2 < handle
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
        LambdaT   %[0, 0, 0, lambda4 - 1/2]
        Sigma0    %Longterm variance
        revlevel  %mean reverting level of the VAR model
        T         %maturities in month
        L4Vec     
        LambdaRN  
        
        % inovs and HistData are passed by reference
        inovs     %inovation terms generate from Normal(0, 1), 
                  %4 dimensions:10000 scenarios for 5 variables for 180 months starting from 301 time points in historical data
        HistData  %Historical Data
        HistRest  %Historical Residual 
        HistHMat  %Historical Ht Matrix
        Sigma     %Sigma from VAR estimation
        HistTStrtData %Historical Data for bond yield at different maturities
        
        boundedMat
    end

    
    methods        
        function obj = RNEstimation_2(suminfo, HistData, Sigma, revlevel,inovs, HistTStrtData)
            % Constructor
            if nargin > 0
                obj.revlevel = revlevel;
                obj.update(suminfo);
                obj.HistData = HistData;
                obj.Sigma = Sigma;
                obj.inovs = inovs;
                obj.T = [12 36 60 84 120 144 168 180];
                obj.HistHMat = suminfo.CondVariances;
                obj.HistRest = suminfo.Residuals;
                obj.L4Vec = [0; 0; 0; suminfo.lambda4];
                %obj.LambdaRN = diag(obj.LambdaT - obj.L4Vec);
                obj.LambdaRN = [0; 0; 0; -0.5];
                obj.HistTStrtData = HistTStrtData(:, obj.T / 12);
                
                obj.boundedMat = gpuArray(25.*repmat(obj.Sigma0,301,10000)); %50
            end
        end

    
        function [mu, Beta, Cmat, A1, B1, LambdaT, Sigma0] = load(obj, suminfo)
            % Load parameter vector into corresponding variables
            Beta = reshape(suminfo.xCenter(1:16),[4,4]);
            Cmat = diag(suminfo.xCenter(17:20));
            A1 = diag(suminfo.xCenter(21:24));
            B1 = diag(suminfo.xCenter(25:28));
            
            LambdaT = [0; 0; 0; suminfo.lambda4-0.5];
            
            mu = (eye(4) - Beta) * obj.revlevel;
            %mu = (eye(5) - Beta) * (obj.revlevel - [-0.43;-0.25;-0.0000;-0.003;-0.027]);
            
            Sigma0 = suminfo.xCenter(17:20)./(ones(4,1)-suminfo.xCenter(21:24)-suminfo.xCenter(25:28));
        end

        function obj = update(obj, suminfo)
            % Update parameters
            [obj.mu, obj.Beta, obj.Cmat, obj.A1, obj.B1, obj.LambdaT, obj.Sigma0] = obj.load(suminfo);
        end
        
        function boundedSigmas = boundedSigma(obj, Sigmas)
                %[~,N] = size(Sigmas);
                boundedSigmas = max(eps,min(Sigmas, obj.boundedMat));
        end
        
        
        function yTermLSE = RNEstGPU(obj, lambda)
            m_idx = (1:4:1204)';
            A = gpuArray(repmat(diag(obj.A1),301,1));
            B = gpuArray(repmat(diag(obj.B1),301,1));
            C = gpuArray(repmat(diag(obj.Cmat),301,1));
            L4V = gpuArray(repmat(obj.L4Vec,301,1));
            LRN = gpuArray(repmat(obj.LambdaRN, 301, 1));
                   
            Z0 = gpuArray(reshape(obj.HistData.m(1:301,:)', 1204, 1));
            Z1 = gpuArray(reshape(obj.HistData.m(2:302,:)', 1204, 1));
            
            lambda0 = [lambda(1); lambda(2); 0; obj.mu(4)];           
            lambda1 = [lambda(3:6);lambda(7:10);0 0 0 0; obj.Beta(4,:)];
          
            L0 = gpuArray(repmat(lambda0, 301, 1));
            L1 = gpuArray(kron(eye(301), lambda1));
            mu_L0 = gpuArray(repmat(obj.mu - lambda0, 301, 1));
            B_L1 = gpuArray(kron(eye(301), obj.Beta - lambda1));
            
            LastResT = gpuArray(reshape(obj.HistRest', 1204, 1));
            LastSigT = gpuArray(reshape(obj.HistHMat', 1204, 1));
            
            SigmaT = gpuArray(C + A .* ((LastResT - L4V .* LastSigT - (L0 + L1 * Z0)).^2)...
                     + B .* LastSigT);
                 
            SigmaT = gpuArray(repmat(SigmaT, 1, 10000));
            LastZT = gpuArray(repmat(Z1, 1, 10000));          
            yStruct = gpuArray(zeros(301, size(obj.T, 2)));

            sumZT = gpuArray(LastZT(m_idx,:));
            idx = 1;
            
            for i=2:180
                Yt = sqrt(SigmaT) .* obj.inovs.m(:,:,i-1);
                Zt =  (SigmaT .* LRN + mu_L0) + B_L1 * LastZT + Yt;           
                sumZT = sumZT + Zt(m_idx,:);
                %disp(sumZT)
                if i == obj.T(idx)
                    yStruct(:,idx) = -log(mean(exp(-sumZT),2))/obj.T(idx)*12;
                    idx = idx + 1;
                end

%                 if sum(sum(isnan(Zt(1,:)))) > 0
%                     disp(i);
%                     yTermLSE = Zt(1:5,:);
%                     return
%                 end
                
                SigmaT = obj.boundedSigma(C + A .* ...
                         ((Yt - SigmaT .* L4V - (L1 * LastZT + L0)).^2) ...
                         + B .* SigmaT);               
                LastZT = Zt;
            end            
            yTermLSE = gather(sum(sum((obj.HistTStrtData - yStruct).^2)));
        end
                
        function outMat = test(obj, t, lambda, ts, scn)
            Z0 = obj.HistData.m(t-1,:)';
            Z1 = obj.HistData.m(t,:)';
            
            lambda0 = [lambda(1); lambda(2); 0; obj.mu(4)];
            lambda1 = [lambda(3:6);lambda(7:10);0 0 0 0; obj.Beta(4,:)];
            
            LastResT = obj.HistRest(t-1,:)';
            LastSigT = obj.HistHMat(t-1,:)';
            
            SigmaT = diag(obj.Cmat) + obj.A1 * ...
                     (((LastResT - diag(obj.L4Vec) * LastSigT - (lambda0 + lambda1 * Z0)).^2)...
                     + obj.B1 * LastSigT);
            SigmaT = repmat(SigmaT, 1, 10000);
            LastZT = repmat(Z1, 1, 10000);    
            
            yStruct = zeros(1, size(obj.T, 2));
            sumZT =LastZT(1,:);            
            ts = ts + 1;
            outMat = zeros(5, ts-1, 2);
            idx = 1;
            for i=2:180
                Yt = sqrt(SigmaT) .* obj.inovs.m(:,:,i-1,t-1);
                
                Zt =  ((SigmaT .* obj.LambdaRN) + (obj.mu - lambda0)) + ...
                      (obj.Beta - lambda1) * LastZT + Yt;
                
                sumZT = sumZT + exp(Zt(1,:));
                
                % keep track of Zt and SigmaT
                outMat(:,i-1,1) = Zt(:,scn);
                outMat(:,i-1,2) = SigmaT(:,scn);
                
                % return as required time step
                if i==ts
                   return
                end
                
                if i == obj.T(idx)
                    yStruct(1,idx) = -log(mean(exp(-sumZT)))/obj.T(idx)*12;
                    idx = idx + 1;
                end
                
                SigmaT = diag(obj.Cmat) + diag(obj.A1) .* ...
                         ((Yt - SigmaT.*obj.L4Vec - (lambda1*LastZT+lambda0)).^2) ...
                         + diag(obj.B1) .* SigmaT;
                
                LastZT = Zt;
            end
        end 
        
        function yTermLSE = testGPU(obj, lambda, ts, scn)
            m_idx = (1:4:1204)';
            A = gpuArray(repmat(diag(obj.A1),301,1));
            B = gpuArray(repmat(diag(obj.B1),301,1));
            C = gpuArray(repmat(diag(obj.Cmat),301,1));
            L4V = gpuArray(repmat(obj.L4Vec,301,1));
            LRN = gpuArray(repmat(obj.LambdaRN, 301, 1));
                   
            Z0 = gpuArray(reshape(obj.HistData.m(1:301,:)', 1204, 1));
            Z1 = gpuArray(reshape(obj.HistData.m(2:302,:)', 1204, 1));
            
            lambda0 = [lambda(1); lambda(2); 0; obj.mu(4)];           
            lambda1 = [lambda(3:6);lambda(7:10);0 0 0 0; obj.Beta(4,:)];
          
            L0 = gpuArray(repmat(lambda0, 301, 1));
            L1 = gpuArray(kron(eye(301), lambda1));
            mu_L0 = gpuArray(repmat(obj.mu - lambda0, 301, 1));
            B_L1 = gpuArray(kron(eye(301), (obj.Beta - lambda1)));
            
            LastResT = gpuArray(reshape(obj.HistRest', 1505, 1));
            LastSigT = gpuArray(reshape(obj.HistHMat', 1505, 1));
            
            SigmaT = gpuArray(C + A .* (((LastResT - L4V .* LastSigT - (L0 + L1 * Z0)).^2)...
                     + B .* LastSigT));
                 
            SigmaT = gpuArray(repmat(SigmaT, 1, 10000));
            LastZT = gpuArray(repmat(Z1, 1, 10000));          
            yStruct = gpuArray(zeros(301, size(obj.T, 2)));
            sumZT = gpuArray(exp(LastZT(m_idx,:)));
            idx = 1;
            
            % Edit later
            %ts = 34;
            ts = ts + 1;
            outMat = zeros(5, ts-1, 3);

            for i=2:180
                Yt = sqrt(SigmaT) .* obj.inovs.m(:,:,i-1);
                temp = B_L1 * LastZT;
                Zt =  (SigmaT .* LRN + mu_L0) + B_L1 * LastZT + Yt;           
                sumZT = sumZT + exp(Zt(m_idx,:));
                if i == obj.T(idx)
                    yStruct(:,idx) = -log(mean(exp(-sumZT),2))/obj.T(idx)*12;
                    idx = idx + 1;
                end
                
                outMat(:,i-1,1) = gather(Zt(1:5,scn));
                outMat(:,i-1,2) = gather(SigmaT(1:5,scn));
                
                if sum(sum(isnan(Zt(1,:)))) > 0
                    disp(i);
                    return
                end
                
                SigmaT = obj.boundedSigma(C + A .* ...
                         ((Yt - SigmaT .* L4V - (L1 * LastZT + L0)).^2) ...
                         + B .* SigmaT);
                
                LastZT = Zt;
            end
            yTermLSE = gather(sum(sum((obj.HistTStrtData - yStruct).^2)));
        end      
        
        function outMat = test2(obj, t, lambda, ts, scn)
            Z0 = obj.HistData.m(t-1,:)';
            Z1 = obj.HistData.m(t,:)';
            
            lambda0 = [lambda(1); lambda(2); 0; obj.mu(4)];
            lambda1 = [lambda(3:6);lambda(7:10);0 0 0 0; obj.Beta(4,:)];
            
            LastResT = obj.HistRest(t-1,:)';
            LastSigT = obj.HistHMat(t-1,:)';
            
            SigmaT = diag(obj.Cmat) + obj.A1 * ...
                     (((LastResT - diag(obj.L4Vec) * LastSigT - (lambda0 + lambda1 * Z0)).^2)...
                     + obj.B1 * LastSigT);
            SigmaT = repmat(SigmaT, 1, 10000);
            LastZT = repmat(Z1, 1, 10000);    
            
            yStruct = zeros(1, size(obj.T, 2));
            sumZT = exp(LastZT(1,:));            
            ts = ts + 1;
            outMat = zeros(5, ts-1, 2);
            idx = 1;
            for i=2:180
                Yt = sqrt(SigmaT) .* obj.inovs.m(:,:,i-1,t-1);
                
                Zt =  ((SigmaT .* obj.LambdaRN) + (obj.mu - lambda0)) + ...
                      (obj.Beta - lambda1) * LastZT + Yt;
                
                sumZT = sumZT + exp(Zt(1,:));
                
                % keep track of Zt and SigmaT
                outMat(:,i-1,1) = Zt(:,scn);
                outMat(:,i-1,2) = SigmaT(:,scn);
                
                % return as required time step
                if i==ts
                   return
                end
                
                if i == obj.T(idx)
                    yStruct(1,idx) = -log(mean(exp(-sumZT)))/obj.T(idx)*12;
                    idx = idx + 1;
                end
                
                SigmaT = obj.boundedSigma( diag(obj.Cmat) + diag(obj.A1) .* ...
                         ((Yt - SigmaT.*obj.L4Vec - (lambda1*LastZT+lambda0)).^2) ...
                         + diag(obj.B1) .* SigmaT );
                
                LastZT = Zt;
            end
        end 
        

        
        %CPU version
        function yStruct = termstrt(obj, t, lambda)
            Z0 = obj.HistData.m(t-1,:)';
            Z1 = obj.HistData.m(t,:)';
            
            lambda0 = [lambda(1); lambda(2); 0; obj.mu(4)];
            lambda1 = [lambda(3:6);lambda(7:10);0 0 0 0; obj.Beta(4,:)];

            LastResT = obj.HistRest(t-1,:)';
            LastSigT = obj.HistHMat(t-1,:)';
            
            SigmaT = diag(obj.Cmat) + diag(obj.A1) .* ...
                     ((LastResT - diag(obj.L4Vec) * LastSigT - (lambda0 + lambda1 * Z0)).^2)...
                     + diag(obj.B1) .* LastSigT;
            SigmaT = repmat(SigmaT, 1, 10000);
            LastZT = repmat(Z1, 1, 10000);     
            
            yStruct = zeros(1, size(obj.T, 2));
            sumZT = LastZT(1,:);
            idx = 1;
            for i=2:180
                Yt = sqrt(SigmaT) .* obj.inovs.m(:,:,i-1,t-1);               
                Zt =  (SigmaT .* obj.LambdaRN) + (obj.mu - lambda0) + ...
                      (obj.Beta - lambda1) * LastZT + Yt;               
                sumZT = sumZT + Zt(1,:);
                
                % check for NaN
%                 if sum(sum(isnan(Zt(1,:)))) > 0
%                     disp(i);
%                     yStruct = Zt;
%                     break
%                 end
                
                if i == obj.T(idx)
                    yStruct(1,idx) = -log(mean(exp(-sumZT)))/obj.T(idx)*12;
                    idx = idx + 1;
                end
                
%                 SigmaT = diag(obj.Cmat) + diag(obj.A1) .* ...
%                          ((Yt - SigmaT.*obj.L4Vec - (lambda1*LastZT+lambda0)).^2) ...
%                          + diag(obj.B1) .* SigmaT;
                  SigmaT = obj.boundedSigma( diag(obj.Cmat) + diag(obj.A1) .* ...
                           ((Yt - SigmaT.*obj.L4Vec - (lambda1*LastZT+lambda0)).^2) ...
                           + diag(obj.B1) .* SigmaT );
                
                  LastZT = Zt;
            end
        end

        function yTermLSE = RNEst(obj, lambda)
            yTermRN = zeros(301, size(obj.T, 2));
            for t=2:302
                yTermRN(t-1,:) = obj.termstrt(t, lambda);
            end
            yTermLSE = sum(sum((obj.HistTStrtData - yTermRN).^2));
            %yTermLSE = sumsqr(obj.HistTStrtData - yTermRN);
        end
        
%           Correct Index AS IN THESIS?
%             sumZT = 0;
%             idx = 1;
%             
%             for i=1:180
%                 Yt = sqrt(SigmaT) .* obj.inovs.m(:,:,i);
%                 Zt =  (SigmaT .* LRN + mu_L0) + B_L1 * LastZT + Yt;           
%                 sumZT = sumZT + exp(Zt(m_idx,:));
%                 if i == obj.T(idx)
%                     yStruct(:,idx) = -log(mean(exp(-sumZT),2))/obj.T(idx)*12;
%                     idx = idx + 1;
%                 end
% 
% %                 if sum(sum(isnan(Zt(1,:)))) > 0
% %                     disp(i);
% %                     yTermLSE = Zt(1:5,:);
% %                     return
% %                 end
%                 
%                 SigmaT = obj.boundedSigma(C + A .* ...
%                          ((Yt - SigmaT .* L4V - (L1 * LastZT + L0)).^2) ...
%                          + B .* SigmaT);               
%                 LastZT = Zt;
%             end
        
    end
end


