clear 
close all
clc
main_matlab_K_1
main_matlab_K_2
%% PLANT MODEL
load identified_t1
load identified_t1_full

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = zeros(1,2);
np = sqrt(numel(Ap));
nu = numel(Bp)/np;
ny = numel(Cp)/np;
save plant Ap Bp Hp Cp Dp

%%
%%%%%%% Kalman Filter Design
sysk=ss(A,[Bp Hp],C,0);

[kest,L,P] = kalman(sysk,5e3,1,0);

%%%%%%% Optimal Gain
syst=ss(A,B,C,0);
[K,S,e] = lqry(syst,1e4,1e0*eye(nu),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ak = rlqg.A;
Bk = rlqg.B;
Ck = rlqg.C;
Dk = rlqg.D;
save controller_reduced.mat Ak Bk Ck Dk

%% PLANT MODEL
load plant
load controller_reduced

nk = sqrt(numel(Ak));

%%
load F1
load F2
load E1
load E2
E = blkdiag(E1,E2);
F = [F1;F2];
[eig(A) eig(A+B*F)]

Mt = ss(A+B*F,B*E^-1,F,E^-1)-E^-1;
Nt = ss(A+B*F,B*E^-1,C+Dp*F,Dp*E^-1);

sim('AW_test')