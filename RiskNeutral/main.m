%% Read Data
%paramRN = csvread("Input/GARCHparam.csv");
%paramRN = csvread("InputXcenter25/xCenter25.csv");

load("..\\Data\\suminfo1");
Sigma = csvread("..\\Data\\Sigma.csv");
HistData = csvread("..\\Data\\HistData.csv");
HistMean  = csvread("..\\Data\\HistMean.csv");
HistTStrtData = csvread("..\\Data\\HistTermStrt.csv");
%lambda = csvread("..\\Data\\lambdaSrDiv25(3).csv");

lambda = [-0.0147770119  0.0104606557, 2.5696250518 -0.0192116997  3.1913911088  0.0076979464 -1.8220007978  0.0007271041 -2.2535170262 -0.0064646368];

%HistRest = csvread("InputXcenter25/HistRest.csv");
%HistHMat  = csvread("InputXcenter25/HistHMat.csv");

% Pass by reference objects
HistData = RefMat(HistData);
inovs = RefMat(randn(1204,10000,179,'single')); %for GPU only
%inovs = RefMat(randn(1505,10000,180,'single'));
%inovs = RefMat(randn(5,10000,179,301,'single'));
%inovs = RefMat(randn(5,10000,661,'single'));

%% Routine for Optimization
RNEstimt = RNEstimation_2(suminfo, HistData, Sigma, HistMean, inovs, HistTStrtData);
options = optimset('Display','iter','MaxIter', 100, 'MaxFunEvals', 100);
%startX = [-0.0012 0.0009 0.0022 0.0006 -0.0015 0.0001 -0.0008 0.0017 0.0003 0.0001 0.0005 -0.0010];
%startX = [-0.00103 0.0009 0.043608 0.0291 0 0.04753 0.004652 -0.01422 0 0.10245];
%startX = [-0.0017 0.0002 0.0057 0.0015 -0.0007 -0.0003 -0.0043 0.0055 -0.0001 0.0002 -0.0000 -0.0041];
startX = [-0.0009 0.0003 0.005068 -0.051777 0.09866 0.24261 0.0174424 -0.04471 -0.052807	0.01047];
%x = fminsearch(@RNEstimt.RNEst, startX, options);
x6 = fminsearch(@RNEstimt.RNEstGPU, startX, options); %%current
startX = x;
x7 = fminsearch(@RNEstimt.RNEstGPU, x6, options);
csvwrite("lambdaSrDiv26(1).csv",x);
x5 = fminsearch(@RNEstimt.RNEstGPU, lambda, options);
x3 = fminsearch(@RNEstimt.RNEstGPU, x_test, options); 
x4 = fminsearch(@RNEstimt.RNEstGPU, x_test, options); 
x2 = fminsearch(@RNEstimt.RNEstGPU, x_test, options); 

%x5 = fminsearch(@RNEstimt.RNEstGPU, lambda, options);
tic
    RNEstimt.RNEstGPU(startX)
toc

%% Routine for Test
RNEstimt = RNEstimation_2(paramRN, HistData, Sigma, HistMean, HistHMat, HistRest, inovs, HistTStrtData);

% will output the first time step that NaN occurs
Outmat = RNEstimt.termstrt(2, zeros(1,12));  
% find the first fail scenario
[row, col] = find(isnan(Outmat));
% total number of fail cases up to the first Inf value
  totalNaN = sum(sum(isnan(Outmat))); 

ts = 60; % display information up to this time step
scnCol = col(1); % pick the fail scenario 

%% Test for CPU version
% get Z_t and Sigma_t for the fail scenario up to ts
outMat = RNEstimt.test(2, zeros(1,12), ts, scnCol);
testZT = outMat(:,:,1);
testSigmaT = outMat(:,:,2);
% get corresponding innovation terms
testInovsT = squeeze(inovs.m(:,scnCol,1:ts,1));

% plot Sigma_t, Z_t, Innovation term
figure
for i=1:5
    subplot(5,3,i*3-2) 
    plot(testSigmaT(i,:));
    
    subplot(5,3,i*3-1) 
    plot(testZT(i,:));
    
    subplot(5,3,i*3)
    plot(testInovsT(i,:));
end

% display information for a normal scenario (i.e. pick scenario 1)
OutmatNormal = RNEstimt.test2(2, zeros(1,12), ts, scnCol);
testZTNormal = OutmatNormal(:,:,1);
testSigmaTNormal = OutmatNormal(:,:,2);
testInovsTNormal = squeeze(inovs.m(:,scnCol,1:ts,1));


figure
for i=1:5
    subplot(5,3,i*3-2) 
    plot(testSigmaTNormal(i,:));
    
    subplot(5,3,i*3-1) 
    plot(testZTNormal(i,:));
    
    subplot(5,3,i*3)
    plot(testInovsTNormal(i,:));
end

%% Test for GPU version
% get Z_t and Sigma_t for the fail scenario up to ts
RNEstimt = RNEstimation_2(paramRN, HistData, Sigma, HistMean, HistHMat, HistRest, inovs, HistTStrtData);
ts = 30;
scnCol = 7067;
outMat = RNEstimt.testGPU(zeros(1,12), ts, scnCol);
B_L1 = [0.980610000000000,0.0430550000000000,1.18700000000000,-0.275620000000000,-0.0191940000000000;0.00543820000000000,0.975570000000000,0.716640000000000,-0.0596360000000000,-0.00221440000000000;-2.11310000000000e-05,-0.000488220000000000,0.124580000000000,0.00190070000000000,-0.00195080000000000;0,0,0,0,0;0.00103690000000000,-0.0128720000000000,-0.153260000000000,-1.02690000000000,0.989350000000000];
lastZT = [-4.96818399429321;-5.00366783142090;-0.00110406638123095;0.0622519329190254;-5.82175970077515];
testZT = outMat(:,:,1);
testSigmaT = outMat(:,:,2);
testInovsT = outMat(:,:,3);
% get corresponding innovation terms
testInovsT = squeeze(inovs.m(1:5,scnCol,1:ts));

% plot Sigma_t, Z_t, Innovation term
figure
for i=1:5
    subplot(5,3,i*3-2) 
    plot(testSigmaT(i,:));
    
    subplot(5,3,i*3-1) 
    plot(testZT(i,:));
    
    subplot(5,3,i*3)
    plot(testInovsT(i,:));
end

% display information for a normal scenario (i.e. pick scenario 1)
normalScn = 1;
OutmatNormal = RNEstimt.test(2, zeros(1,12), ts, normalScn);
testZTNormal = OutmatNormal(:,:,1);
testSigmaTNormal = OutmatNormal(:,:,2);
testInovsTNormal = squeeze(inovs.m(:,normalScn,1:ts,1));

figure
for i=1:5
    subplot(5,3,i*3-2) 
    plot(testSigmaTNormal(i,:));
    
    subplot(5,3,i*3-1) 
    plot(testZTNormal(i,:));
    
    subplot(5,3,i*3)
    plot(testInovsTNormal(i,:));
end


%% Show Final Result
outMat = RNEstimt.RNEstGPU(x2);
outMat = gather(outMat);
actMat = RNEstimt.HistTStrtData;
MSEmat = (outMat - actMat).^2;
totalError = sum(sum(MSEmat));
MisPricTerm = outMat-actMat;
mean_MisPric = mean(MisPricTerm,1);
std_MisPric = var(MisPricTerm,1).^0.5;

opt1Lambda0 = [x2(1);x2(2);0;RNEstimt.mu(4)];
opt1Lambda1 = [x2(3:6); x2(7:10); 0 0 0 0 ; RNEstimt.Beta(4,:)];
RNBeta = RNEstimt.Beta-opt1Lambda1;

plotX = (1:1:301);
figure
plot(plotX, actMat(:,8),plotX, outMat(:,8));

plotM = datetime(1991,5,1) + calmonths(0:300);
plotT = [12 36 60 84 120 144 168 180];
h = figure;
for i = 1:8
    subplot(4,2,i)
    plot(plotM, actMat(:,i),'k-', plotM, outMat(:,i), 'r-.')
    %xlim(datetime( 2016,[5 6],[1 1]))
    %xtickformat('dd-MMM-yyyy')
    legend('Historical Yield', 'Fitted Yield')
    title(strcat("Annual Yield of ", convertCharsToStrings(num2str(plotT(i)/12)),"-Year Bond")); %
    xlabel("Year");
end

sim_result = RNEstimt.generator(HistMean, x2);

plotM = datetime(2016,6,1) + calmonths(0:660);
names = {'Short-Rate','Long-Rate','Inflation','Stock Returns'};
fig6 = figure('units','pixels','outerposition',[0 0 1000 800]);
for di = 1:4
    var1 = squeeze(sim_result.Zt(:,di,:))';
    prct1 = prctile(var1,[5 25 50 75 95], 1);
    subplot(2,2,di)
    plot(plotM, prct1(1,2:662), plotM, prct1(2,2:662), plotM, prct1(3,2:662), plotM, prct1(4,2:662), plotM, prct1(5,2:662));
    title(names{di});
end
suptitle('Funnels of Doubt');

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'rnfit','-dpdf','-r0')



trow = floor(301/12);
RMSEtable = zeros(trow, 8);
for i = 1:trow
    RMSEtable(i,:) = sum(MSEmat((i*12-11):i*12,:).^0.5);
end
RMSEtable(trow,:) = RMSEtable(trow,:) + MSEmat(301, :).^0.5;
