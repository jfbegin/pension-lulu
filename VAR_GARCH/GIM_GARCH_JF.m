classdef GIM_GARCH_JF < handle
  % GIM_GARCH
  
  properties
    series = [];
    HistMean = NaN(4,1);      % Average level of the series (?)
    Sigma = NaN(4,1);         % Variance level with classic VAR model
    
    % VAR parameters
    Mu      = NaN(4,1);
    Beta    = NaN(4,4);
    LambdaT = NaN(4,1);
    
    % GARCH parameters
    Omega  = NaN(4,4);
    Alpha1 = NaN(4,4);
    Beta1  = NaN(4,4);
    
    xCenter = NaN
    
    % !!! These should be private, i.e., should not be changed.
    numData = NaN;
    Sigma0  = NaN;
  end
  
  methods
    function obj = GIM_GARCH_JF(param, series, Sigma, HistMean)
    % Constructor
      if nargin == 4
        obj.HistMean = HistMean;
        obj.series = series;
        obj.numData = size(obj.series,1);
        obj.Sigma = Sigma;
        obj.update(param);
      else 
        error('Bad number of input parameters!');
      end
    end
    
    function [Mu, Beta, Omega, Alpha1, Beta1, LambdaT, Sigma0] = load(obj, param)
    % Load param into their respective variables
      Beta   = reshape(param(1:16),[4,4]);
      Mu     = (eye(4) - Beta) * (obj.HistMean);
      
      Omega  = diag(param(17:20));
      Alpha1 = diag(param(21:24));
      Beta1  = diag(param(25:28));
      
      Sigma0 = Omega./(1-Alpha1-Beta1);
      LambdaT = [0; 0; 0; mean(obj.series(:,4))/Sigma0(4,4)];
    end
    
    function obj = update(obj, param)
    % Update parameters
      [obj.Mu, obj.Beta, obj.Omega, obj.Alpha1, obj.Beta1, obj.LambdaT, obj.Sigma0] = obj.load(param);
    end
    
    function logliket = logLLH_sub(obj, res, h)
      logliket = -2 * log(2 * pi) - 0.5 * log(det(h)) - 0.5 * res * pinv(h) * res';
      % logliket = log(mvnpdf(res,zeros(1,4),h));
    end
    
    function boundedSigmas = boundedSigma(obj, h)
      boundedSigmas = max(eps, min(h, 25.*obj.Sigma0));
    end    
      
    function boundedSigmas = boundedSigmaScen(obj, Sigmas)
      %[~,N] = size(Sigmas);
      boundedMatScen = 25.*repmat(diag(obj.Sigma0), 1, 10000);
      boundedSigmas = max(eps,min(Sigmas, boundedMatScen));
    end
    
    function [logL,Resi,H] = negLogLLHCenter(obj, param)
    % Function for calculating log likelihood
      
      % Update parameters
      update(obj, param);

      if max(obj.Alpha1 + obj.Beta1) > 0.99
      	logL = -Inf;
      	return
      end

      if min([diag(obj.Omega); diag(obj.Alpha1); diag(obj.Beta1)]) < 1e-10
        logL = -Inf;
        return
      end
      
      % Check stationarity for the VAR model
      if max(abs(eig(obj.Beta))) > 0.9957
        logL = -Inf;
        return
      end
      
      H(1,:) = diag(obj.Sigma0);
      Ht = diag(H(1,:));
      Resi(1,:) = obj.series(2,:)' - obj.Beta * obj.series(1,:)' - Ht * obj.LambdaT;
      
      logL(1) = obj.logLLH_sub(Resi(1,:), Ht);
      
      for dt = 3:obj.numData
        Et = obj.Beta * obj.series(dt - 1,:)';
        
        resmatterm = diag(diag(obj.Alpha1) .* (Resi(dt-2,:).^2));
        htmatterm = diag(diag(obj.Beta1) .* (H(dt-2,:)));
        H(dt-1,:) = diag(obj.boundedSigma(diag(obj.Omega + diag(resmatterm + htmatterm))));

        Ht = diag(H(dt-1,:));
        Resi(dt-1,:) = obj.series(dt,:)' - Et - Ht * obj.LambdaT;

        logL(dt-1) =  obj.logLLH_sub(Resi(dt-1,:),Ht);
      end
    end
    
    function suminfo = optimize(obj, paraminit)
    % Routine that optimizes the model
      optionsset = optimset('Display','iter','MaxIter',100000, 'MaxFunEvals', 100000);
      obj.xCenter = fminsearch(@(x) -sum(obj.negLogLLHCenter(x)), paraminit, optionsset);

      [logL,Resi,H] = obj.negLogLLHCenter(obj.xCenter);
      update(obj, obj.xCenter);
      
      suminfo.Residuals = Resi;
      suminfo.CondVariances = H;
      suminfo.LL = sum(logL);
    end
    
    function sim_result = generator(obj, startz0)
      %Simulate real-world scenarios
      
      sim_result.Zt = NaN(663,4,10000);
      sim_result.Resi = NaN(662,4,10000);
      sim_result.Ht = NaN(662,4,10000);
      
      sim_result.Zt(1,:,:) = repmat(startz0, 1, 10000);
      
      Z1 = startz0;
            
      rng('default');
      inovt = randn(4,10000,661,'single');
            
%     LastResT = obj.HistRest(301,:)';
%     LastSigT = obj.HistHMat(301,:)';
            
      LastResi = diag(obj.Sigma0).^(1/2);
      LastHt = diag(obj.Sigma0);
                     
      Sigmat = diag(obj.Omega) + diag(obj.Alpha1).*(LastResi.^2) + diag(obj.Beta1) .* LastHt;           
      Sigmat = repmat(Sigmat, 1, 10000);
      
      sim_result.Ht(1,:,:) = repmat(LastHt, 1, 10000);  
      sim_result.Zt(2,:,:) = repmat(Z1, 1, 10000); 
      sim_result.Resi(1,:,:) = repmat(LastResi, 1, 10000);
      
      for i=2:662
        %disp(i);
        sim_result.Ht(i,:,:) = Sigmat;
        Yt = sqrt(Sigmat) .* inovt(:,:,i-1);   
        sim_result.Resi(i,:,:) = Yt;
                
        %disp(size(diag(obj.Beta) .* obj.simZt(i-1,:,:)));
        %disp(size(Yt));
        %disp(size(obj.LambdaT .* SigmaT));
        Zt =  obj.Mu + obj.Beta* squeeze(sim_result.Zt(i,:,:)) + Yt + obj.LambdaT .* Sigmat; % modified

        Sigmat = obj.boundedSigmaScen(diag(obj.Omega) + diag(obj.Alpha1) .* (Yt.^2) + diag(obj.Beta1) .* Sigmat);       
                
        sim_result.Zt(i+1,:,:) = Zt;
      end
    end

  end
end

