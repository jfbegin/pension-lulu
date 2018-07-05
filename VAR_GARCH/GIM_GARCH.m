classdef GIM_GARCH < handle
    % GIM_GARCH
    
    properties
        paramVec
        mu
        Beta
        Cmat
        A1
        B1
        LambdaT
        numData
        Sigma0
        revlevel
        
        lb
        ub
        % pointer to data, accessed by GARCHinput.m
        GARCHinput_
        Sigma
    end
    
    methods        
        function obj = GIM_GARCH(param, GARCHinput, Sigma, revlevel)
            % Constructor
            if nargin > 0
                obj.revlevel = revlevel;
                obj.update(param);
                obj.GARCHinput_ = GARCHinput;
                obj.numData = size(GARCHinput.m, 1);
                obj.Sigma = Sigma;                
            end
        end
        
        function [mu, Beta, Cmat, A1, B1, LambdaT, Sigma0] = load(obj, param)
            % Load parameter vector into corresponding variables
            Beta = reshape(param(1:16),[4,4]);
            %disp(obj.Beta);
            %disp(obj.revlevel);
            mu = (eye(4) - Beta) * (obj.revlevel);
            
            Cmat = diag(param(17:20));
            A1 = diag(param(21:24));
            B1 = diag(param(25:28));

            Sigma0 = param(17:20)./(ones(4,1)-param(21:24)-param(25:28));
            LambdaT = [0; 0; 0; 0.003893229/Sigma0(4)-0.5];

        end

        function obj = update(obj, param)
            % Update parameters
            [obj.mu, obj.Beta, obj.Cmat, obj.A1, obj.B1, obj.LambdaT, obj.Sigma0] = obj.load(param);
        end
        
        function logliket = logLLH_sub(~, ht, residualt)
            % Helper function for calculating log likelihood
            %disp(ht);
            %disp(residualt);
            logliket = -2 * log(2 * pi) - 0.5 * log(det(ht)) - 0.5 * residualt' * pinv(ht) * residualt;
        end       

        
        function boundedSigmas = boundedSigma(~, Sigmas, Sigma0)
            %[~,N] = size(Sigmas);
             boundedSigmas = max(eps,min(Sigmas, 25.*Sigma0));
        end
        
        function loglikeT = negLogLLHCenter(obj, param)
            % Function for calculating log likelihood
            %[mu, Beta, Cmat, A1, B1, LambdaT] = obj.load(param);
            
            if min(param(17:28)) < 10 ^ (-11)
                %disp("garch positive");
                %disp(min(param(26:41)));
                loglikeT = Inf;
                return
            end

            
            for i=(1:4)
               if param(20+i) + param(24+i) >= 0.99
                   %disp("garch stationary");
                   %disp(param(30+i) + param(35+i));
                   %disp(i);
                   loglikeT = Inf;
                   return
               end
            end
            
            %load parameters
            Beta = reshape(param(1:16),[4,4]);
            Cmat = diag(param(17:20));
            A1 = diag(param(21:24));
            B1 = diag(param(25:28));
            Sigma0 = param(17:20)./(ones(4,1)-param(21:24)-param(25:28));
 %           LambdaT = [-0.5; -0.5; 0; param(41)-0.5; -0.5];
            LambdaT = [0; 0; 0; 0.003893229/Sigma0(4)-0.5];   
            
            %check stationarity for the VAR model
            if max(abs(eig(Beta))) >=0.9957
                %disp("VAR stationary");
                %disp(max(abs(eig(Beta))));
                loglikeT = Inf;
                return
            end
            
            %load data
            GARCHinput = obj.GARCHinput_;
            %Htm2 = diag(diag(obj.Sigma));
            
            Htm2 = diag(Sigma0);
            Residualtm2 = GARCHinput.m(2,:)' - Beta * GARCHinput.m(1,:)' - Htm2 * LambdaT;
            
            loglikeT = obj.logLLH_sub(Htm2, Residualtm2);
            
            for t = 3:obj.numData
                Et = Beta * GARCHinput.m(t - 1,:)';

                resmatterm = diag(diag(A1) .* (Residualtm2.^2));
                htmatterm = diag(diag(B1) .* diag(Htm2));
                %Ht = Cmat + resmatterm + htmatterm;
                
                %Try to put a bound there
                Ht = Cmat + resmatterm + htmatterm;
                Ht = diag(obj.boundedSigma(diag(Ht),Sigma0));
                
                %disp(t);
                %disp(GARCHinput.m(t,:));
                Residualt = GARCHinput.m(t,:)' - Et - Ht * LambdaT;
                
%                 disp(t);
%                 disp(Ht);
%                 disp(Residualt);
%                 disp(pinv(Ht));
%                 disp('------------');
                loglikeT = loglikeT + obj.logLLH_sub(Ht, Residualt);
                Residualtm2 = Residualt;
                Htm2 = Ht;
            end
            loglikeT = -loglikeT;
        end
        
    end
end

