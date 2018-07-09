%% Read Data
paramRN = csvread("InputXcenter25/xCenter25.csv");
Sigma = csvread("InputXcenter25/Sigma.csv");
HistData = csvread("InputXcenter25/HistData.csv");
HistMean = csvread("InputXcenter25/HistMean.csv");
HistRest = csvread("InputXcenter25/HistRest.csv");
HistHMat  = csvread("InputXcenter25/HistHMat.csv");
HistTStrtData = csvread("InputXcenter25/HistTermStrt.csv");
lambda = csvread("lambdaSrDiv25(3).csv");

% Pass by reference objects
HistData = RefMat(HistData);
%inovs = RefMat(randn(5,10000,179,301,'single'));
%inovs = RefMat(randn(5,10000,661,'single'));
%clearvars inovs
rng('default');
%inovs = RefMat(randn(3305,10000,179,'single')); %for GPU only
inovs = RefMat(randn(280,10000,179,'single'));
%inovs = RefMat(randn(3305,10000,179,'single')); %for GPU only
%inovs = RefMat(randn(5,10,179,301,'single'));

%% Routine for Optimization
rnScen = SimResultRN(paramRN, HistData, Sigma, HistMean, HistHMat, HistRest, inovs, HistTStrtData, lambda);
%outMat = rwScen.RNEstGPU(lambda);
rnScen.genScenario();
outMat = rnScen.oneSceBondYieldYear(1);
rnScen.fillTermStruct(10000);
%outMat = gather(outMat);

rn_simHt= rnScen.simHt;
rn_simRest = rnScen.simRest;
rn_simZt = rnScen.simZt;
rn_termStruct = rnScen.termStruct;

save('rn_simHt.mat', 'rn_simHt');
save('rn_simRest.mat', 'rn_simRest');
save('rn_simZt.mat', 'rn_simZt');
save('rn_termStruct.mat', 'rn_termStruct');

[X,Y,Z] = ind2sub(size(rnScen.termStruct(:,:,1:200)), find(rnScen.termStruct(:,:,1:200)==0));
test_scen = unique(Z);
[X,Y,Z] = ind2sub(size(rnScen.termStruct(:,:,1:200)), find(isinf(rnScen.termStruct(:,:,1:200))));
test2_scen = unique(Z);

plotX = (1:1:662);
figure;

for i = 1:5
    var1 = squeeze(rnScen.simZt(:,i,:))';
    prct1 = prctile(var1,[5 25 50 75 95], 1);
    subplot(3,2,i)
    plot(plotX, prct1(1,2:663), plotX, prct1(2,2:663), plotX, prct1(3,2:663), plotX, prct1(4,2:663), plotX, prct1(5,2:663));
end

figure;

plot(plotX, squeeze(outMat(:,1)))

figure;

for i = 1:5
    var1 = squeeze(rnScen.simHt(:,i,200))';
    subplot(3,2,i);
    plot(plotX, var1);
end

figure;
for i = 1:5
    var1 = squeeze(rnScen.simZt(2:663,i,200))';
    subplot(3,2,i);
    plot(plotX, var1);
end

