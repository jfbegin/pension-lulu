load("../Data/rn_simHt.mat"); %rn_simHt = temp;
load("../Data/rn_simRest.mat"); %rn_simRest = temp;
load("../Data/rn_simZt.mat"); %rn_simZt = temp;
load("../Data/rn_termStruct.mat"); %rn_termStruct = outMat;

load("../Data/rw_simHt.mat"); %rn_simHt = temp;
load("../Data/rw_simRest.mat"); %rn_simRest = temp;
load("../Data/rw_simZt.mat"); %rn_simZt = temp;
load("../Data/rw_termStruct.mat"); %rn_termStruct = outMat;


%Get the 0 and inf scenarios for both risk neutral and real world
[X,Y,Z] = ind2sub(size(rn_termStruct), find(rn_termStruct==0));
rn0Scen = unique(Z);

[X,Y,Z] = ind2sub(size(rw_termStruct), find(rw_termStruct==0));
rw0Scen = unique(Z);

[X,Y,Z] = ind2sub(size(rn_termStruct), find(isinf(rn_termStruct)));
rnInfScen = unique(Z);

[X,Y,Z] = ind2sub(size(rw_termStruct), find(isinf(rw_termStruct)));
rwInfScen = unique(Z);

%united all the bad scenarios
infScen = union([rnInfScen;rwInfScen],[rn0Scen;rw0Scen]);
%infScen = union(rwInfScen, rw0Scen);

%clean out the bad scenarios
goodScen = setdiff(1:10000, infScen);


% [X,Y,Z] = ind2sub(size(rw_termStruct), find(isinf(rw_termStruct)));
% nanScen = unique(Z);
% goodScen = setdiff([1:1:10000], nanScen);

rw_termStruct = rw_termStruct(:,:,goodScen);
rw_simHt = rw_simHt(:,:,goodScen);
rw_simRest = rw_simRest(:,:,goodScen);
rw_simZt = rw_simZt(:,:,goodScen);

rw_Zt = rw_simZt;

rn_termStruct = rn_termStruct(:,:,goodScen);
rn_simHt = rn_simHt(:,:,goodScen);
rn_simRest = rn_simRest(:,:,goodScen);
rn_simZt = rn_simZt(:,:,goodScen);

rn_Zt = rn_simZt;

%%Plot of Projected Scenarios
plotM = datetime(2016,6,1) + calmonths(0:660);
h = figure;

plotTitle = ["Short-term Bond Yield", "Long-term Bond Yield", "Inflation Rate", "Excess Stock Return", "Dividend Yield"];

for i = 1:4
    var1 = squeeze(rw_Zt(:,i,:))';
    prct1 = prctile(var1,[5 25 50 75 95], 1);
    subplot(3,2,i)    
    plot(plotM, prct1(1,2:662), 'k-.',  plotM, prct1(2,2:662), 'k--', plotM, prct1(3,2:662), 'k-',  plotM, prct1(4,2:662), 'k--', plotM, prct1(5,2:662), 'k-.');
    title(plotTitle(i));
    xlabel("Year");
    ylabel("Monthly Rate")
    legend('5-th Percentile','25-th Percentile','50-th Percentile','75-th Percentile','95-th Percentile');   
end

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'projrate','-dpdf','-r0')

%% Plot of Historical Data
% Dep. of plot calculation
if false
  plotM = datetime(1991,5,1) + calmonths(0:301);
  h = figure;

  plotHist(:,[1 2 5]) = exp(HistData(:,[1 2 5]));

  plotTitle = ["Short-term Bond Yield", "Long-term Bond Yield", "Inflation Rate", "Excess Stock Return", "Dividend Yield"];

  for i = 1:5
      subplot(5,1,i)    
      plot(plot( :,i));
      title(plotTitle(i));
      xlabel("Year");
      ylabel("Monthly Rate")
  end

  set(h,'Units','Inches');
  pos = get(h,'Position');
  set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
  print(h,'projrate','-dpdf','-r0')
  
  plotX = (1:1:661);
  figure;

  for i = 1:4
      var1 = squeeze(rn_Zt(:,i,:))';
      prct1 = prctile(var1,[5 25 50 75 95], 1);
      subplot(3,2,i)
      plot(plotX, prct1(1,2:662), plotX, prct1(2,2:662), plotX, prct1(3,2:662), plotX, prct1(4,2:662), plotX, prct1(5,2:662));
  end
end

%% Get the Asset Scenarios
rn_assetScenario = AssetScenario(rn_Zt, size(goodScen,2), 55, rn_termStruct);
rw_assetScenario = AssetScenario(rw_Zt, size(goodScen,2), 55, rw_termStruct);

%% Clean Asset Return > 200
% call ALMStudy first -- rn_ALMStudy = ALMStudy(55, 9306, rn_assetScenario, 1, TargetCalAssump, PlanDesign);

%% STOP !!! We need to run the liability code before running this part. 

a = rn_ALMStudy.ARMat;
[X,Y] = ind2sub(size(a), find(a>200));
TooBigIdx = unique(Y);

b = rw_ALMStudy.ARMat;
[X,Y] = ind2sub(size(b), find(b>200));
TooBigIdxrw = unique(Y);
%TooBigIdx = csvread('TooBigAR.csv');
%goodScen_new = setdiff(1:size(goodScen_bound,2),TooBigIdx);
goodScen_new = setdiff(1:size(goodScen,2),TooBigIdx);

rw_assetScenario2 = AssetScenario(rw_Zt(:,:,goodScen_new), size(goodScen_new,2), 55, rw_termStruct(:,:,goodScen_new));
rn_assetScenario2 = AssetScenario(rn_Zt(:,:,goodScen_new), size(goodScen_new,2), 55, rn_termStruct(:,:,goodScen_new));

%% Bounded Mat Cleaning
%assetScenario = AssetScenario(n_Zt, 9986, 55, n_termStruct);
% goodScen_bound = goodScen;
% 
% Sigma0 = param(26:30)./(ones(5,1)-param(31:35)-param(36:40));
% boundedMat = 25.*repmat(Sigma0', 661, 1, size(goodScen_bound,2));
% 
% test = rw_simHt > boundedMat + 0.00001;
% 
% [X,Y,Z] = ind2sub(size(test), find(test));
% nanScen = unique(Z);
% goodScen_bound = setdiff(1:1:size(goodScen_bound,2), nanScen);
% 
% rw_termStruct = rw_termStruct(:,:,goodScen_bound);
% rw_simHt = rw_simHt(:,:,goodScen_bound);
% rw_simRest = rw_simRest(:,:,goodScen_bound);
% rw_simZt = rw_simZt(:,:,goodScen_bound);
% 
% rw_Zt = rw_simZt;
% rw_Zt(:,[1 2 5],:) = exp(rw_simZt(:,[1 2 5],:));
% 
% rn_termStruct = rn_termStruct(:,:,goodScen_bound);
% rn_simHt = rn_simHt(:,:,goodScen_bound);
% rn_simRest = rn_simRest(:,:,goodScen_bound);
% rn_simZt = rn_simZt(:,:,goodScen_bound);
% 
% rn_Zt = rn_simZt;
% rn_Zt(:,[1 2 5],:) = exp(rn_simZt(:,[1 2 5],:));
% 
% 
% %% After clean all the bad scens, rerun the the asset Scenario
% rw_assetScenario = AssetScenario(rw_Zt, size(goodScen_bound,2), 55, rw_termStruct);
% rn_assetScenario = AssetScenario(rn_Zt, size(goodScen_bound,2), 55, rn_termStruct);



% sum(sum(isinf(rw_assetScenario.srAnnualEfft)))
% sum(sum(isinf(rw_assetScenario.divAnnualEfft)))
% sum(sum(isinf(rw_assetScenario.bondAnnualEfft)))
% sum(sum(isinf(rw_assetScenario.infAnnualEfft)))
% 
% max(max(rw_assetScenario.bondYield_S_Efft))
% max(max(rw_assetScenario.bondYield_L_Efft))
% 
% 
% sum(sum(isinf(rn_assetScenario.srAnnualEfft)))
% sum(sum(isinf(rn_assetScenario.divAnnualEfft)))
% sum(sum(isinf(rn_assetScenario.bondAnnualEfft)))
% sum(sum(isinf(rn_assetScenario.infAnnualEfft)))
% 
% max(max(rn_assetScenario.bondYield_S_Efft))
% max(max(rn_assetScenario.bondYield_L_Efft))