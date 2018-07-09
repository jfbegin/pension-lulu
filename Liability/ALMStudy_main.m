%% Plan Design Code

TargetCalAssump = struct('ValRateType', 'BondBased', 'StockWeight',0.5, 'BondWeight', 0.5, 'EquityRP', 0.0223, ...
                         'Lbound', 0.01, 'Ubound', 0.1, 'FixValRate', 0.05, 'SetContri', false, 'TargetBen', 0.01, 'ContriRate', 0.2);
                                          

PlanDesign = struct('ValMethod','PUC','ValRateType','BondBased','EquityRP',0.023,'StockWeight',0.5,'BondWeight',0.5,...
                    'FixValRate',0.05,'Lbound',0.01,'Ubound',0.1,'NoActionRange',false,'TriggerLFR',1,'TriggerUFR',1.4,...
                    'AdjLFR',1,'AdjUFR',1.4,'AdjAction','edge');
        
PlanDesignEAN = struct('ValMethod','PUC','ValRateType','EROA','EquityRP',0.023,'StockWeight',0.5,'BondWeight',0.5,...
                    'FixValRate',0.05,'Lbound',0.01,'Ubound',0.1,'NoActionRange',true,'TriggerLFR',1,'TriggerUFR',1.4,...
                    'AdjLFR',1,'AdjUFR',1.4,'AdjAction','edge', 'ActionFP','FP');
                
PlanDesignDC = struct('ValRateType','BondBased','EquityRP',0.023,'StockWeight',0.5,'BondWeight',0.5,...
                    'FixValRate',0.05,'Lbound',0.01,'Ubound',0.1);  
                
                
%% ALM Study Code                
% rw_ALMStudy = ALMStudy(55, 9976, rw_assetScenario, 1, TargetCalAssump, PlanDesign);
% rn_ALMStudy = ALMStudy(55, 9976, rn_assetScenario, 1, TargetCalAssump, PlanDesign);
% 
% rn_DCStudy = DCALMStudy(55,9976,rn_assetScenario,TargetCalAssump,PlanDesignDC);
% rw_DCStudy = DCALMStudy(55,9976,rw_assetScenario,TargetCalAssump,PlanDesignDC);
% 
% rn_VBALMStudy = VBALMStudy(rn_ALMStudy,rn_DCStudy);

rw_ALMStudy = ALMStudyEAN(55, 9976, rw_assetScenario, 1, TargetCalAssump, PlanDesignEAN);
rn_ALMStudy = ALMStudyEAN(55, 9976, rn_assetScenario, 1, TargetCalAssump, PlanDesignEAN);

rn_DCStudy = DCALMStudy(55,9976,rn_assetScenario,TargetCalAssump,PlanDesignDC);
rw_DCStudy = DCALMStudy(55,9976,rw_assetScenario,TargetCalAssump,PlanDesignDC);

rn_VBALMStudy = VBALMStudyEAN(rn_ALMStudy,rn_DCStudy);



%% Result AND Check Consistency
rw_ALMStudy.InitialLiaIC
rw_DCStudy.InitialFund

%rw_ALMStudy.InitialTargetPUC
%rw_ALMStudy.InitialLiaPUC

rw_BondReturn = rw_assetScenario.bondAnnualEfft;
rw_StockReturn = rw_assetScenario.srAnnualEfft;
rw_TermStruct = rw_assetScenario.termStruct;

rw_TermStrucn1Scen = rw_TermStruct(:,:,1);

TBPAR = rw_ALMStudy.ARMat;
DCAR = rw_DCStudy.ARMat;

TBPDiscount=rw_ALMStudy.ValRateMat;
DCDiscount = rw_DCStudy.ValRateMat;

TBPBenefit = rw_ALMStudy.SimulationResult.AdjAccBen;
DCBenefit = rw_DCStudy.SimulationResult.BenPaid;

TBPBenefit1Scen = TBPBenefit(:,:,1);
DCBenefit1Scen = DCBenefit(:,:,1);

TBPSalary = rw_ALMStudy.salary;
DCSalary = rw_DCStudy.salary;

TBPReplaceRatio = rw_ALMStudy.ReplaceRatioAtRetire;
DCReplaceRatio = rw_DCStudy.ReplaceRatioAtRetire;

TBPFundRatio = rw_ALMStudy.SimulationResult.FundRatio;

TBPSalary1Scen = TBPSalary(:,:,1);
DCSalary1Scen = DCSalary(:,:,1);

TBPGenAccount = rn_VBALMStudy.CashflowBasic;
DCGenAccount = rn_VBALMStudy.DCCashflowBasic;

TBPGenAccount1Scen = TBPGenAccount(:,:,1);
DCGenAccount1Scen = DCGenAccount(:,:,1);

TBPAccrualRate = rw_ALMStudy.SimulationResult.AccrualRate;
TBPAccrualRate31 = squeeze(rw_ALMStudy.SimulationResult.AccrualRate(2,:,:));
TBPAccrualRate64 = squeeze(rw_ALMStudy.SimulationResult.AccrualRate(35,:,:));

%Generational Value 
TBPGenValue = rn_VBALMStudy.VGenAccountBasic - mean(rn_VBALMStudy.InitialContriPaid, 2);
DCGenValue = rn_VBALMStudy.VDCGenAccountBasic - mean(rn_VBALMStudy.DCInitialContriPaid, 2);

%Surplus Deficit Option
GenSurplus = rn_VBALMStudy.GenSurplus;
GenDeficit = rn_VBALMStudy.GenDeficit;
GenResSurplus = rn_VBALMStudy.GenResSurplus;
GenResDeficit = rn_VBALMStudy.GenResDeficit;
GenTotalSurplus = rn_VBALMStudy.GenTotalSurplus;
GenTotalDeficit = rn_VBALMStudy.GenTotalDeficit;

%Value Transfer
Vtransfer = rn_VBALMStudy.VGenAccountBasic - BaseResult.VGenAccountBasic;
BVtransfer = rn_VBALMStudy.VGenAccountBen - BaseResult.VGenAccountBen;
RVtransfer = rn_VBALMStudy.VGenAccountRes - BaseResult.VGenAccountRes;



%% Plot
plotX55 = (1:1:(rw_ALMStudy.ProjYear));
plotX56 = (1:1:(rw_ALMStudy.ProjYear+1));
plotGen = (-25:85);

h=figure;
plot(plotX55, squeeze(TBPAR(:,1)), '-.', plotX55, squeeze(TBPAR(:,3)), '-');
legend('Path 1', 'Path 2')
ylabel("Fund Return")
xlabel("Years")
ylim([-0.2,0.35]);
xlim([0,55]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'fr1scen','-dpdf','-r0')


h = figure;
prctAR = prctile(TBPAR',[5 25 50 75 95], 1);
plot(plotX55, prctAR(1,:), 'k-.', plotX55, prctAR(2,:),'k--', plotX55, prctAR(3,:),'k-', plotX55, prctAR(4,:),'k--', plotX55, prctAR(5,:),'k-.');
legend('5-th Percentile','25-th Percentile','50-th Percentile','75-th Percentile','95-th Percentile')
ylabel("Fund Return")
xlabel("Years");
ylim([-0.2,0.3]);
xlim([0,55]);


set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'frpct1','-dpdf','-r0')


%Plot of key operation variables : valuation rate, spread between valuation
%rate and fund return, fund ratio, projected accrual rate

h = figure;
%subplot(3,1,1);
prctDCT = prctile(TBPDiscount',[5 25 50 75 95], 1);
plot(plotX56, prctDCT(1,:), 'k-.', plotX56, prctDCT(2,:),'k--', plotX56, prctDCT(3,:),'k-', plotX56, prctDCT(4,:),'k--', plotX56, prctDCT(5,:),'k-.');
ylabel("Valuation Rate")
%title("Distribution of valuation rate")
xlabel("Years");
legend('5-th Percentile','25-th Percentile','50-th Percentile','75-th Percentile','95-th Percentile')
xlim([0,55]);
ylim([0.015,0.09]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'DISC1','-dpdf','-r0')


%subplot(3,1,2);
% AVSpread = TBPAR' - TBPDiscount(1:55,:)';
% prctAVS = prctile(AVSpread,[5 25 50 75 95], 1);
h = figure;
prctAVS = prctAR - prctDCT(:,1:55);
plot( plotX55, prctAVS(2,:),'k--', plotX55, prctAVS(3,:),'k-', plotX55, prctAVS(4,:),'k--');
ylabel("Spread")
%legend('25-th percentile', 'median', '75-th percentile')
%title("Distribution of spread between investment return and valuation rate")
xlabel("Years");
legend('25-th Percentile','50-th Percentile','75-th Percentile')
xlim([0,55]);
ylim([-0.085,0.08]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Spread1','-dpdf','-r0')

subplot(3,1,3);
h = figure;
prctFR = prctile(TBPFundRatio',[5 25 50 75 95], 1);
plot(plotX56, prctFR(1,:), 'k-.', plotX56, prctFR(2,:),'k--', plotX56, prctFR(3,:),'k-', plotX56, prctFR(4,:),'k--', plotX56, prctFR(5,:),'k-.');
legend('5-th Percentile','25-th Percentile','50-th Percentile','75-th Percentile','95-th Percentile')
%title("Distribution of projected funded ratio")
ylabel("Funded Ratio")
xlabel("Years");
xlim([0 55])
ylim([0.7 1.8])

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'FR3','-dpdf','-r0')

% Plot of Accrual Rate Dynamics
h=figure;
plot(plotX56, prctile(squeeze(TBPAccrualRate(2,:,:))',50, 1), plotX56, prctile(squeeze(TBPAccrualRate(12,:,:))',50, 1), plotX56, prctile(squeeze(TBPAccrualRate(22,:,:))',50, 1), plotX56, prctile(squeeze(TBPAccrualRate(32,:,:))',50, 1),plotX56, prctile(squeeze(TBPAccrualRate(42,:,:))',50, 1));
legend('age 31', 'age 41', 'age 51', 'age 61', 'age 71')
title("Median accrual rate at different age")
ylabel("Accrual rate")
xlabel("Years");

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'MedianAccR3','-dpdf','-r0')


% Plot of Retirement Ratio
%Plot of TBP Replacement Ratio at Retirement
h=figure;
prctRR = prctile(TBPReplaceRatio,[5 25 50 75 95], 1);
plot(plotX55, prctRR(1,:), 'k-.', plotX55, prctRR(2,:),'k--', plotX55, prctRR(3,:),'k-', plotX55, prctRR(4,:),'k--', plotX55, prctRR(5,:),'k-.', plotX55, ones(55).*0.010 * 35, 'r-');

ylabel("Replacement Ratio")
xlabel("Years");
legend('5-th Percentile','25-th Percentile','50-th Percentile','75-th Percentile','95-th Percentile', 'Initial Target')
xlim([0,55]);
ylim([0.05,1]);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'RR','-dpdf','-r0')

%Plot of DC Replacement Ratio at Retirment
h = figure;

prctDCRR = prctile(DCReplaceRatio,[5 25 50 75 95], 1);
plot(plotX55, prctDCRR(1,:), 'k-.', plotX55, prctDCRR(2,:),'k--', plotX55, prctDCRR(3,:),'k-', plotX55, prctDCRR(4,:),'k--', plotX55, prctDCRR(5,:),'k-.', plotX55, ones(55).*0.0103 * 35, 'r-');

ylabel("Replacement Ratio")
xlabel("Years");


%Plot of Benefit Adjustment Distribution
%figure                                                                                
%bar(rw_ALMStudy.DistReturnList.ProbIncDec','stack');
h = figure;
%bar(rw_ALMStudy.DistReturnList.ProbIncDec','stack');
colorb = bar(rw_ALMStudy.DistReturnList.DistBenAdj','stack', 'FaceColor','flat');
for k = 1:9
   colorb(k).CData = k;
end

ylabel("Distribution of Benefit Adjustments")
xlabel("Years");
ylim([0 1]);
xlim([0 55]);

%barl = legend('10%+', '5-10%+', '2-5%+', '0-2%+', 'no change', '0-2%-', '2-5%-', '5-10%-', '10%-');
labels = {'10%-', '5-10%-', '2-5%-', '0-2%-', 'no change', '0-2%+', '2-5%+', '5-10%+', '10%+'};
fliplegend(labels);

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'ProbIncDc3','-dpdf','-r0')

%'b-', plotGen, BaseResult.VGenAccountBasic - mean(BaseResult.InitialContriPaid, 2), 

%Plot Value of Pension Deal
h = figure;
plot(plotGen, TBPGenValue, 'b-', plotGen, BaseResult.VGenAccountBasic - mean(BaseResult.InitialContriPaid, 2), 'k-.', plotGen, ones(111).*0, 'r-');
xlim([-25 85]);
ylabel("Value of the Pension Deal")
xlabel("Cohorts");
legend("Saving Corridor Plan", "EROA Plan")

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'VPension3','-dpdf','-r0')

%Plot of Option Values
h = figure;
subplot(3,1,1);
plot(plotGen,GenSurplus, 'b-', plotGen,GenDeficit, 'r-.');
xlim([-25 85]);
ylim([0 65000]);
ylabel("Option Value")
xlabel("Cohorts");
legend("VBC", "VBP")
title("Benefit Options")

subplot(3,1,2);
plot(plotGen,GenResSurplus, 'b-', plotGen, GenResDeficit, 'r-.');
xlim([-25 85]);
ylim([0 65000]);
ylabel("Option Value")
xlabel("Cohorts");
legend("VRC", "VRP")
title("Residual Options")

subplot(3,1,3);
plot(plotGen,GenTotalSurplus,'b-', plotGen,GenTotalDeficit,'r-.');
xlim([-25 85]);
ylim([0 65000]);
ylabel("Option Value")
xlabel("Cohorts");
legend("Call Options", "Put Options")
title("Total Value of Put and Call Options")

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'SDOptions1','-dpdf','-r0')

%Plot of V-transfer
h=figure;
subplot(2,1,1);
plot(plotGen, Vtransfer, 'k-', plotGen, ones(111).*0, 'r-')
xlim([-25 85])
title("Total Value Transfer")
ylabel("Value Transfer")
xlabel("Cohorts");

subplot(2,1,2);
plot(plotGen, BVtransfer, 'b-.', plotGen, RVtransfer, 'k--',plotGen, ones(111).*0, 'r-')
xlim([-25 85])
legend('Benefit Value', 'Residual Value')
title("Benefit and Residual Value Transfer")
ylabel("Value Transfer")
xlabel("Cohorts");

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Vtrans3','-dpdf','-r0')

%% Other Plot
%Plot of stock return and bond return percentile
figure
prctBR = prctile(rw_BondReturn',[5 25 50 75 95], 1);
plot(plotX55, prctBR(1,:), plotX55, prctBR(2,:), plotX55, prctBR(3,:), plotX55, prctBR(4,:), plotX55, prctBR(5,:))

figure
prctSR = prctile(rw_StockReturn',[5 25 50 75 95], 1);
plot(plotX55, prctSR(1,:), plotX55, prctSR(2,:), plotX55, prctSR(3,:), plotX55, prctSR(4,:), plotX55, prctSR(5,:));

%Plot of stock return and bond return for one scenario
figure
plot(plotX55, squeeze(rw_BondReturn(:,1)), '-.', plotX55,  squeeze(rw_StockReturn(:,1)),'-');

%Plot of Term Struct 1 Scen
figure
plot((1:1:15), rw_TermStrucn1Scen(3,:), (1:1:15), rw_TermStrucn1Scen(24,:),(1:1:15), rw_TermStrucn1Scen(50,:));

%Plot of First Year Ret Benefit
figure

prctFirstB = prctile(squeeze(TBPBenefit(36,:,:))',[5 25 50 75 95], 1);
plot(plotX56, prctFirstB(1,:), plotX56, prctFirstB(2,:), plotX56, prctFirstB(3,:), plotX56, prctFirstB(4,:), plotX56, prctFirstB(5,:));

%Plot of Last Year Salary
figure
prctLastS = prctile(squeeze(TBPSalary(64,:,:))',[5 25 50 75 95], 1);
plot(plotX55, prctLastS(1,:), plotX55, prctLastS(2,:), plotX55, prctLastS(3,:), plotX55, prctLastS(4,:), plotX55, prctLastS(5,:));


%Plot risk neutral Fund Ratio Distribution
figure
prctrnFR = prctile(rn_ALMStudy.SimulationResult.FundRatio',[5 25 50 75 95], 1);
plot(plotX56, prctrnFR(1,:), plotX56, prctrnFR(2,:), plotX56, prctrnFR(3,:), plotX56, prctrnFR(4,:), plotX56, prctrnFR(5,:));

figure
c = rn_VBALMStudy.VDCGenAccountBasic - mean(rn_VBALMStudy.DCInitialContriPaid, 2);
plot(plotGen, c);

%% EROA
% TargetCalAssumpS = struct('ValRateType', 'EROA', 'StockWeight', 0.5, 'BondWeight', 0.5, 'EquityRP', 0.02, ...
%                          'Lbound', 0.01, 'Ubound', 0.1, 'FixValRate', 0.05, 'SetContri', false, 'TargetBen', 0.01, 'ContriRate', 0.2);
% 
% PlanDesignS = struct('ValMethod','PUC','ValRateType','EROA','EquityRP',0.02,'StockWeight',0.5,'BondWeight',0.5,...
%                     'FixValRate',0.05,'Lbound',0.01,'Ubound',0.1,'NoActionRange',true,'TriggerLFR',0.9,'TriggerUFR',1.1,...
%                     'AdjLFR',0.9,'AdjUFR',1.1,'AdjAction','edge');
% PlanDesignS1 = struct('ValMethod','PUC','ValRateType','EROA','EquityRP',0.02,'StockWeight',0.5,'BondWeight',0.5,...
%                     'FixValRate',0.05,'Lbound',0.01,'Ubound',0.1,'NoActionRange',false,'TriggerLFR',0.9,'TriggerUFR',1.1,...
%                     'AdjLFR',0.9,'AdjUFR',1.1,'AdjAction','edge');
% PlanDesignDCS = struct('ValRateType','bondbased','EquityRP',0.02,'StockWeight',0.5,'BondWeight',0.5,...
%                     'FixValRate',0.05,'Lbound',0.01,'Ubound',0.1);                 
%                 
% rn_ALMStudyS = ALMStudy(55, 9031, rn_assetScenario2, 1, TargetCalAssumpS, PlanDesignS);
% rn_DCStudyS = DCALMStudy(55,9031,rn_assetScenario2,TargetCalAssumpS,PlanDesignDC);
% rn_ALMStudyS1 = ALMStudy(55, 9031, rn_assetScenario2, 1, TargetCalAssumpS, PlanDesignS1);
% rn_VBALMStudyS = VBALMStudy(rn_ALMStudyS,rn_DCStudyS);
% rn_VBALMStudyS1 = VBALMStudy(rn_ALMStudyS1,rn_DCStudyS);
% 
% %a = vb_studyResult.GenAccountScen;
% aS = rn_VBALMStudyS.VGenAccountBasic - mean(rn_VBALMStudyS.InitialContriPaid, 2);
% Vtransfer = rn_VBALMStudyS1.VGenAccountBasic - rn_VBALMStudyS.VGenAccountBasic;
% 
% cS = rn_VBALMStudyS.VDCGenAccountBasic - mean(rn_VBALMStudyS.DCInitialContriPaid, 2);
% 
% figure;
% plot(-24:85, cS)
% 
% hold on
% plot(1, 100000,'w');
% 
% 
% figure;
% plot(-24:85, aS);
% 
% GenSurplusS = rn_VBALMStudy.GenSurplus;
% GenDeficitS = rn_VBALMStudy.GenDeficit;
% GenResSurplusS = rn_VBALMStudy.GenResSurplus;
% GenResDeficitS = rn_VBALMStudy.GenResDeficit;
% GenTotalSurplusS = rn_VBALMStudy.GenTotalSurplus;
% GenTotalDeficitS = rn_VBALMStudy.GenTotalDeficit;
% 
% 
% 
% plotX = (-24:85);
% figure;
% subplot(3,2,1);
% plot(plotX,GenSurplusS);
% subplot(3,2,2);
% plot(plotX,GenDeficitS);
% subplot(3,2,3);
% plot(plotX,GenResSurplusS);
% subplot(3,2,4);
% plot(plotX,GenResDeficitS);
% subplot(3,2,5);
% plot(plotX,GenTotalSurplusS);
% subplot(3,2,6);
% plot(plotX,GenTotalDeficitS);
% 
 BaseResult = rn_VBALMStudy;