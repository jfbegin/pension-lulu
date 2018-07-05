
%% Main Routain

% load initialization values
%paramCenter = csvread("D:\\Fall Research\\Matlab\\VAR_GARCH\\xCenter25.csv");
paramCenter = csvread("..\\Data\\Starting_param.csv");

% load historical data
GARCHinputCen = csvread("..\\Data\\HistData.csv");
GARCHinputCen = RefMat(GARCHinputCen);

% VAR sigma
Sigma = csvread("..\\Data\\Sigma.csv");

% historical mean
HistMean = csvread("..\\Data\\HistMean.csv");


%param_test = param+rand(45,1);
%obj = GIM_GARCH(param, GARCHinput, Sigma);
objCenter =  GIM_GARCH(paramCenter, GARCHinputCen, Sigma, HistMean);

options2 = optimset('Display','iter','MaxIter',100000, 'MaxFunEvals', 100000);
%options3 = optimset('Display','iter','MaxIter',100000, 'MaxFunEvals', 100000, 'TolX', 0.0000000000000000000000000001);

xCenter = fminsearch(@objCenter.negLogLLHCenter, paramCenter, options2);

%[xCenter2,xfval,xexitflag,xoutput,xgrad,xhessian] = fminunc(@objCenter.negLogLLHCenter, testparam, options3);


param29 = 0.003893/(xCenter(20)/(1-xCenter(24)-xCenter(28)));
xCenter1 = [xCenter; param29];
csvwrite("xCenter.csv",xCenter1);

csvwrite("xCenter27.csv",xCenter1);
csvwrite("xCenter7.csv",paramMlh);

(xCenter1(29))*(xCenter1(20)/(1-xCenter(24)-xCenter(28)));


%% Ignore this part!
%Different Starting Value Trial

%loop for trying different starting values
mlh = 10;
paramMlh = paramCenter;

for i = 1:100
    paramCenter(1:25) = param(5:29);
    paramCenter(20) = -0.9;
    paramCenter(25) = 0.94;
    
    paramCenter(26:30) = diag(Sigma)/50;
    paramCenter(41) = 1.5;
    
    for j = 1:5
        paramCenter(30+j) = 0.1+(0.5-0.1)*rand;
        paramCenter(35+j) = 0.3 + ((0.99 - paramCenter(30+j))-0.3)*rand;
    end
    
    xCenter = fminsearch(@objCenter.negLogLLHCenter, paramCenter, options2);
    currMlh = objCenter.negLogLLHCenter(xCenter);
    
    if currMlh<mlh 
        mlh = currMlh;
        paramMlh = xCenter;
    end
end

    
    

x2 = fminsearch(@obj.negLogLLH, x2, options2);

x3 = fminsearch(@obj.negLogLLH, param, options2);

%% Stationary check
for i=1:5
    disp(x2(34+i) + x2(39+i));
%     if x2(34+i) + x2(39+i) >= 0.995
%         disp(i);
%     end
end


%% Starting value setting
% try GARCH(1,1) for initialization
param_test = param;
Mdl = arima('ARLags',1,'Variance',garch(1,1));
for i=1:5
    EstMdl = estimate(Mdl,GARCHinput.m(:,i));
    param_test(29+i) = EstMdl.Variance.Constant;
    param_test(34+i) = EstMdl.Variance.ARCH{1};
    param_test(39+i) = EstMdl.Variance.GARCH{1};
end

MdlCenter = arima('Constant',0,'ARLags',1,'Variance',garch(1,1));
for i=1:5
    EstMdlCenter = estimate(MdlCenter,GARCHinputCen.m(:,i));
    paramCenter(25+i) = EstMdlCenter.Variance.Constant;
    paramCenter(30+i) = EstMdlCenter.Variance.ARCH{1};
    paramCenter(35+i) = EstMdlCenter.Variance.GARCH{1};
end
paramCenter(26) = 0.0002;
paramCenter(27) = 0.0001;
paramCenter(31) = 0.3;
paramCenter(32) = 0.3;
paramCenter(33) = 0.15;
paramCenter(36) = 0.6;
paramCenter(37) = 0.65;
paramCenter(39) = 0.6;
paramCenter(38) =0.8;

paramCenter(41) = 1.5;
 

paramCenter = xCenter;
paramCenter(20) = -0.9;
paramCenter(25) = 0.94;


    