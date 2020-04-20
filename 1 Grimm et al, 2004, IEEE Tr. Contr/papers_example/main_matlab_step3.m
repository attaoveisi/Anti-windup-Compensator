clear 
close all
clc

%% Load variables
load variables_for_step3
load plant
load controller

nc = nn;
nu = 1;

%% LMI Definition
setlmis([]);

GAMA = lmivar(2,[naw+nv naw+nu]);

%LMI terms
lmiterm([1 1 1 GAMA],H',G,'s'); % LMI #1:
lmiterm([1 1 1 0],Sai); % LMI #1:

LMISYS = getlmis;

%% Finding Feasible Solution 
% [~,x] = feasp(LMISYS,[0 1000 -1 99 0]);
[~,x] = feasp(LMISYS);
GAMA = dec2mat(LMISYS,x,1);
GAMA1 = GAMA(1:naw,1:naw);
GAMA2 = GAMA(1:naw,naw+1:end);
GAMA3 = GAMA(naw+1:end,1:naw);
GAMA4 = GAMA(naw+1:end,naw+1:end);

save GAMA GAMA1 GAMA2 GAMA3 GAMA4