load('D:\\Lulu\\rn_termStruct.mat');
load('D:\\Lulu\\rn_simZt.mat');
rnZt = temp;
rnTermStruct = outMat;

load('D:\\Lulu\\simZt.mat');
load('D:\\Lulu\\termStruct_total.mat');
rwZt = temp;
rwTermStruct = outMat;

[X,Y,Z] = ind2sub(size(rnTermStruct), find(rnTermStruct==0));
rn0Scen = unique(Z);

[X,Y,Z] = ind2sub(size(rwTermStruct), find(rwTermStruct==0));
rw0Scen = unique(Z);

[X,Y,Z] = ind2sub(size(rnTermStruct), find(isinf(rnTermStruct)));
rnInfScen = unique(Z);

[X,Y,Z] = ind2sub(size(rnZt), find(isinf(rnZt)));
rnZtInfScen = unique(Z);

[X,Y,Z] = ind2sub(size(rwTermStruct), find(isinf(rwTermStruct)));
rwInfScen = unique(Z);

[X,Y,Z] = ind2sub(size(rwZt), find(isinf(rwZt)));
rwZtInfScen = unique(Z);

infScen = union([rnInfScen;rwInfScen],[rn0Scen;rw0Scen;rnZtInfScen;rwZtInfScen]);
goodScen = setdiff(1:10000, infScen);
