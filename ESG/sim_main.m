%% Read Data
%paramRN = csvread("Input/GARCHparam.csv");
load("..\\Data\\suminfo1");
load("..\\Data\\x4_better")

Sigma = csvread("..\\Data\\Sigma.csv");
HistData = csvread("..\\Data\\HistData.csv");
HistMean  = csvread("..\\Data\\HistMean.csv");
HistTStrtData = csvread("..\\Data\\HistTermStrt.csv");

% Pass by reference objects
HistData = RefMat(HistData);
%inovs = RefMat(randn(5,10000,179,301,'single'));
%inovs = RefMat(randn(5,10000,661,'single'));
%clearvars inovs

rng('default');
%inovs = RefMat(randn(3305,10000,179,'single')); %for GPU only
inovs = RefMat(randn(224,10000,179,'single')); %for GPU only
%inovs = RefMat(randn(5,10,179,301,'single'));

%% Routine for Optimization
rwScen = SimResult(suminfo, HistData, Sigma, HistMean, inovs, HistTStrtData, x4);
%outMat = rwScen.RNEstGPU(lambda);
rwScen.genScenario();

% small test
% test = rwScen.simHt > 25.*repmat(rwScen.Sigma0', 662, 1, 10000)+0.00001;
% [X,Y,Z] = ind2sub(size(test), find(test));
% nanStep = unique(X);

% test = squeeze(rwScen.simHt(10,:,:)) > rwScen.boundedMatScen;
% [X,Y] = ind2sub(size(test), find(test));
% nanStep = unique(X);

%outMat = rwScen.oneSceBondYield(1);
%outMat = rwScen.oneSceBondYieldYear(1);
rwScen.fillTermStruct(10000);
%outMat = gather(outMat);

test_Zt = rwScen.simZt;
meanZ1 = mean(squeeze(test_Zt(:,1,:)),2);
meanZ2 = mean(squeeze(test_Zt(:,2,:)),2);
meanZ3 = mean(squeeze(test_Zt(:,3,:)),2);
meanZ4 = mean(squeeze(test_Zt(:,4,:)),2);

lastmeanrate = [meanZ1(663); meanZ2(663); meanZ3(663); meanZ4(663);meanZ5(662)];
max(max(squeeze(test_Zt(:,1,:))))

plotM = datetime(2016,6,1) + calmonths(0:660);


fig6 = figure('units','pixels','outerposition',[0 0 1000 800]);
names = {'Short-Rate','Long-Rate','Inflation','Stock Returns'};
for di = 1:4
    var1 = squeeze(test_Zt(:,di,:))';
    prct1 = prctile(var1,[5 25 50 75 95], 1);
    subplot(2,2,di)
    plot(plotM, prct1(1,2:662), plotM, prct1(2,2:662), plotM, prct1(3,2:662), plotM, prct1(4,2:662), plotM, prct1(5,2:662));
    title(names{di});
end
suptitle('Funnels of Doubt');



figure;

plot(plotX, squeeze(outMat(:,1)))

figure;

for i = 1:5
    var1 = squeeze(rwScen.simHt(:,i,7217))';
    subplot(3,2,i);
    plot(plotX, var1);
end

figure;
for i = 1:5
    var1 = squeeze(rwScen.simZt(2:662,i,7217))';
    subplot(3,2,i);
    plot(plotX, var1);
end

%[X,Y,Z] = ind2sub(size(outMat), find(outMat==0));
%test_scen = unique(Z);
%[X,Y,Z] = ind2sub(size(outMat), find(isinf(outMat)));