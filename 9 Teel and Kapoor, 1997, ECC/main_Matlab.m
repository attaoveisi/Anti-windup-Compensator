clear 
close all
clc

%% PLANT MODEL
% load identified
% load identified_freq
% load identified_t
% B = B/1e4;
load identified_t1
load identified_t1_full

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = zeros(1,2);
save plant Ap Bp Hp Cp Dp

np = sqrt(numel(Ap));
nu = numel(Bp)/np;
nw = numel(Hp)/np;
ny = numel(Cp)/np;

%% Controller
%%%%%%% Kalman Filter Design
sysk = ss(Ap,[Bp Hp],Cp,0);

[kest,L,P] = kalman(sysk,5e3,1,0);

%%%%%%% Optimal Gain
syst = ss(Ap,Bp,Cp,0);
[K,S,e] = lqry(syst,1e4,1e0*eye(nu),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ak = rlqg.A;
Bk = rlqg.B;
Ck = rlqg.C;
Dk = rlqg.D;
save controller.mat Ak Bk Ck Dk
nk = sqrt(numel(Ak));

%%
Q = 1e3*Cp'*Cp;
P = lyap(Ap,Q);

