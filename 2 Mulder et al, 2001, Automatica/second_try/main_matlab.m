clear 
close all
clc

%% PLANT MODEL
load identified_t1
load identified_t1_full

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = zeros(1,2);
save plant Ap Bp Hp Cp Dp

%%
%%%%%%% Kalman Filter Design
sysk=ss(A,[Bp Hp],C,0);

[kest,L,P] = kalman(sysk,5e3,1,0);

%%%%%%% Optimal Gain
syst=ss(A,B,C,0);
[K,S,e] = lqry(syst,1e4,1e0*eye(2),0);

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

np = sqrt(numel(Ap));
nk = sqrt(numel(Ak));
nu = numel(Bp)/np;

%%
Bpu = Bp;
Cpy = Cp;
Cpz = Cpy;
Bpw = Hp;

np = sqrt(numel(Ap));
nk = sqrt(numel(Ak));
nc = nk;
nu = numel(Bpu)/np;
nw = numel(Bpw)/np;
ny = numel(Cpy)/np;
nz = numel(Cpz)/np;

Dpyu = zeros(ny,nu);
Dpyw = zeros(ny,nw);
Dpzu = zeros(nz,nu);
Dpzw = zeros(nz,nw);

plant = ss(Ap,[Bpu,Bpu,Bpw],[Cpz;Cpy],[Dpzu,Dpzu,Dpzw;Dpyu,Dpyu,Dpyw]);
controller = ss(Ak,[Bk,[eye(nk) zeros(nk,nu)]],Ck,[Dk,[zeros(nu,nk) eye(nu)]]);

systemnames = 'plant controller';
inputvar = '[v{2};w;ksi{14}]';
outputvar = '[controller;plant(1)]';
input_to_plant = '[v;controller;w]';
input_to_controller = '[plant(2);ksi]';
cleanupsysic = 'yes';
sys_cl = sysic;
                                                                                                                                                                      
A = sys_cl.A;
nA = sqrt(numel(A));
Bv = sys_cl.B(:,1:nu);
Bw = sys_cl.B(:,nu+1:nu+nw);
Bksi = sys_cl.B(:,nu+nw+1:end);
Cu = sys_cl.C(1:nu,:);
Cz = sys_cl.C(nu+1:end,:);
Duv = sys_cl.D(1:nu,1:nu);
Duw = sys_cl.D(1:nu,nu+1:nu+nw);
Duksi = sys_cl.D(1:nu,nu+nw+1:end);
Dzv = sys_cl.D(nu+1:end,1:nu);
Dzw = sys_cl.D(nu+1:end,nu+1:nu+nw);
Dzksi = sys_cl.D(nu+1,nu+nw+1:end);

%% LMI Definition
setlmis([]);

Q = lmivar(1,[nA 1]);
M = lmivar(1,[1 0;1 0]);
X = lmivar(2,[np+nu nu]);
delta1 = lmivar(2,[1 1]);
Gama = lmivar(1,[1 1]);
% gama = lmivar(1,[1 1]);
gama = 1; % !!!! ATTENTION: gama is not defined as always here. Here Gama is the variable. 

%LMI terms
lmiterm([1 1 1 Q],A,1,'s'); % LMI #1: AQ+Q*A'
lmiterm([1 1 2 0],Bw); % LMI #1: Bw
lmiterm([1 1 3 M],Bv,1); % LMI #1: Bv*M
lmiterm([1 1 3 X],-Bksi,1); % LMI #1: -Bksi*X
lmiterm([1 1 3 Q],1,Cu'); % LMI #1: Q*Cu'
lmiterm([1 1 4 Q],1,Cz'); % LMI #1: Q*Cz'
lmiterm([1 2 2 Gama],-1,1); % LMI #1: -Gama
lmiterm([1 2 3 0],Duw'); % LMI #1: Duw'
lmiterm([1 2 4 0],Dzw'); % LMI #1: Dzw'
lmiterm([1 3 3 M],-2,1); % LMI #1: -2M
lmiterm([1 3 3 M],Duv,1,'s'); % LMI #1: Duv*M+M*Duv'
lmiterm([1 3 3 X],-Duksi,1,'s'); % LMI #1: Duksi*X+X'*Duksi'
lmiterm([1 3 4 M],1,Dzv'); % LMI #1: M*Dzv'
lmiterm([1 3 4 -X],-1,Dzksi'); % LMI #1: -X'*Dzksi'
lmiterm([1 3 5 M],1,1); % LMI #1: M
lmiterm([1 4 4 0],-gama); % LMI #1: -gama*I
lmiterm([1 5 5 delta1],-1,1); % LMI #1: -delta1*I

lmiterm([-2 1 1 Q],1,1); % LMI #3: Q>0

lmiterm([-3 1 1 delta1],1,1); % LMI #4: delta1>0

lmiterm([-4 1 1 Gama],1,1); % LMI #4: Gama>0
% lmiterm([-4 1 1 0],-2); % LMI #4: Gama>0

lmiterm([-5 1 1 M],1,1); % LMI #4: M>0

% lmiterm([-6 1 1 gama],1,1); % LMI #4: M>0

lmiterm([7 1 1 Gama],1,1); % LMI #4: M>0
lmiterm([7 1 1 0],-1.5); % LMI #4: M>0

lmiterm([6 1 1 delta1],1,1); % LMI #4: M>0
lmiterm([6 1 1 0],-1e0); % LMI #4: M>0

LMISYS = getlmis;

%% Finding Feasible Solution 
% Q,M,X,delta1,Gama,gama
n_dec = decnbr(LMISYS);
c = zeros(1,n_dec);
c(end) = 1;
% c(end-1) = 1;
options = [1e-9 1000 0 0 0];
% [~,x] = mincx(LMISYS,c,options,[],[]);
% [~,x] = feasp(LMISYS,[0 1000 -1 99 0]);
[~,x] = feasp(LMISYS);
Q = dec2mat(LMISYS,x,1);
M = dec2mat(LMISYS,x,2);
X = dec2mat(LMISYS,x,3);
delta1 = dec2mat(LMISYS,x,4);
Gama = dec2mat(LMISYS,x,5);
% gama = dec2mat(LMISYS,x,6)

Lambda = X/M
Khat = ss(Ak,[Bk eye(nk) zeros(nk,nu)],Ck,[Dk zeros(nu,nk) eye(nu)]);
sim('LQG_AWBT')