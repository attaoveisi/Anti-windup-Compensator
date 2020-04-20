clear 
close all
clc

%% PLANT MODEL
load variables_for_step2_prt2

%%
% ncl = np+nc;
% nvars = naw*ncl;
% lb = -ones(1,nvars)*1e21;
% ub = ones(1,nvars)*1e20;
% PopulationSize_Data = 10;
% MigrationInterval_Data = 5;
% MigrationFraction_Data = 0.55;
% TolFun_Data = 1e-8;
% PopInitRange_Data = [-1e10;1e10];
% CrossoverFraction_Data = 0.9;
% Generations_Data = 1000;
% [x,fval,exitflag,output,population,score] = ga_opt(nvars,lb,ub,PopInitRange_Data,PopulationSize_Data,CrossoverFraction_Data,MigrationInterval_Data,MigrationFraction_Data,Generations_Data,TolFun_Data);

%% Other method
INS = inv(S);
[a,b] = qr((R*inv(S)*R-R))