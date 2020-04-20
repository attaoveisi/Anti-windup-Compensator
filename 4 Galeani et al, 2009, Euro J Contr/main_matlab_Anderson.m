clear 
close all
clc

%%
load identified_t1

Ap = A;
Bpu = B;
Bpw = H;
Cp = C;
Dpu = zeros(1,2);
Dpw = zeros(1,1);
Cz = C;
Dzu = Dpu;
Dzw = zeros(1,1);

np = sqrt(numel(Ap));
m = numel(Bpu)/np;
q = numel(Bpw)/np;
p = numel(Cp)/np;
l = p;

save plant Ap Bpu Bpw Cp Dpu Dpw Cp Dzw Dzu

%% Controller
%%%%%%% Kalman Filter Design
sysk = ss(Ap,[Bpu Bpw],Cp,0);

[kest,L,P] = kalman(sysk,1e3,1,0);

%%%%%%% Optimal Gain
syst = ss(Ap,Bpu,Cp,0);
[K,S,e] = lqry(syst,1e4,1e0*eye(m),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ac = rlqg.A;
Bc = rlqg.B;
Cc = rlqg.C;
Dc = rlqg.D;
Bcw = zeros(np,q);
Dcw = zeros(m,q);
nc = sqrt(numel(Ac));
save controller.mat Ac Bc Cc Dc Bcw Dcw