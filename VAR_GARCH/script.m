clear all;
clc;

paraminit = csvread("..\\Data\\Starting_param.csv");
series    = csvread("..\\Data\\HistData.csv");
Sigma     = csvread("..\\Data\\Sigma.csv");
HistMean  = csvread("..\\Data\\HistMean.csv");

model =  GIM_GARCH_JF(paraminit, series, Sigma, HistMean);

% suminfo = model.optimize(paraminit);
load('suminfo');
model.update(suminfo.xCenter);
%suminfo.lambda4 = model.LambdaT(4)+1/2;

sim_result = model.generator(HistMean);

plotM = datetime(2016,6,1) + calmonths(0:660);

names = {'Short-Rate','Long-Rate','Inflation','Stock Returns'};
fig1 = figure('units','pixels','outerposition',[0 0 1000 800]);
for di = 1:4
  subplot(2,2,di)
  plot(1:(model.numData),model.series(:,di),'k-');
  xlim([1,model.numData]);
  title(names{di});
end
suptitle('Centred Series');

fig2 = figure('units','pixels','outerposition',[0 0 1000 800]);
for di = 1:4
  subplot(2,2,di)
  plot(1:(model.numData-1),suminfo.Residuals(:,di),'k-');
  xlim([1,model.numData-1]);
  title(names{di});
end
suptitle('Residuals');

fig3 = figure('units','pixels','outerposition',[0 0 1000 800]);
for di = 1:4
  subplot(2,2,di)
  plot(1:(model.numData-1),sqrt(suminfo.CondVariances(:,di)),'k-');
  xlim([1,model.numData-1]);
  title(names{di});
end
suptitle('Conditional Volatilities');

fig4 = figure('units','pixels','outerposition',[0 0 1000 800]);
for di = 1:4
  subplot(2,2,di)
  plot(1:(model.numData-1),suminfo.Residuals(:,di)./sqrt(suminfo.CondVariances(:,di)),'k-');
  xlim([1,model.numData-1]);
  title(names{di});
end
suptitle('Standardized Residuals');

fig5 = figure('units','pixels','outerposition',[0 0 1000 800]);
for di = 1:4
  subplot(2,2,di)
  qqplot(suminfo.Residuals(:,di)./sqrt(suminfo.CondVariances(:,di)));
  title(names{di});
  box on
end
suptitle('QQ-plots');

fig6 = figure('units','pixels','outerposition',[0 0 1000 800]);
for di = 1:4
    var1 = squeeze(sim_result.Zt(:,di,:))';
    prct1 = prctile(var1,[5 25 50 75 95], 1);
    subplot(2,2,di)
    plot(plotM, prct1(1,2:662), plotM, prct1(2,2:662), plotM, prct1(3,2:662), plotM, prct1(4,2:662), plotM, prct1(5,2:662));
    title(names{di});
end
suptitle('Funnels of Doubt');




