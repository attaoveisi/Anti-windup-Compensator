clear 
close all
clc

%% PLANT MODEL
load identified

%%
A = A(1:4,1:4);
B = B(1:4,:);
C = C(:,1:4);
H = H(1:4,:);

Bw=H;
sys1=ss(A,Bw,C,0);

%%%%%%% Kalma Filter Design
sysk=ss(A,[B Bw],C,0);

[kest,L,P] = kalman(sysk,1e3,1,0);

%%%%%%% Optimal Gain
syst=ss(A,B,C,0);
[K,S,e] = lqry(syst,1e5,1e0*eye(2),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ak = rlqg.A;
Bk = rlqg.B;
Ck = rlqg.C;
Dk = rlqg.D;
save controller.mat Ak Bk Ck Dk

%% PLANT MODEL
load plant
load controller

Ap = Ap(1:4,1:4);
Bpu = Bp(1:4,:);
Cpy = Cp(:,1:4);
Cpz = Cpy;
Bpw = Hp(1:4,:);

np = sqrt(numel(Ap));
nk = sqrt(numel(Ak));
nu = numel(Bpu)/np;
nw = numel(Bpw)/np;
ny = numel(Cpy)/np;
nz = numel(Cpz)/np;

Dpyu = zeros(ny,nu);
Dpyw = zeros(ny,nw);
Dpzu = zeros(ny,nu);
Dpzw = zeros(ny,nw);

plant = ss(Ap,[Bpw,Bpu],[Cpz;Cpy],[Dpzw,Dpzu;Dpyw,Dpyu]);
controller = ss(Ak,[zeros(np,nw),Bk],Ck,[zeros(nu,nw),Dk]);

systemnames = 'plant controller';
inputvar = '[w]';
outputvar = '[plant(1)]';
input_to_plant = '[w;controller]';
input_to_controller = '[w;plant(2)]';
% cleanupsysic = 'yes';
sys_cl = sysic;

Acl = sys_cl.A;
Bclw = sys_cl.B;
Cclz = sys_cl.C;
Dclzw = sys_cl.D;

%% LMI Definition
setlmis([]);

R11 = lmivar(1,[np 1]);
R12 = lmivar(2,[np np]);
R22 = lmivar(1,[np 1]);
S = lmivar(1,[np+np 1]);
Gama = lmivar(1,[1 1]);

%LMI terms
lmiterm([1 1 1 R11],1,Ap','s'); 
lmiterm([1 1 2 0],Bpw); 
lmiterm([1 1 3 R11],1,Cpz'); 
lmiterm([1 2 2 Gama],-1,1); 
lmiterm([1 2 3 0],Dpzw'); 
lmiterm([1 3 3 Gama],-1,1); 

lmiterm([2 1 1 S],1,Acl','s'); 
lmiterm([2 1 2 0],Bclw); 
lmiterm([2 1 3 S],1,Cclz'); 
lmiterm([2 2 2 Gama],-1,1); 
lmiterm([2 2 3 0],Dclzw'); 
lmiterm([2 3 3 Gama],-1,1);

lmiterm([-3 1 1 R11],1,1); 
lmiterm([-3 1 2 R12],1,1); 
lmiterm([-3 2 2 R22],1,1); 

lmiterm([-4 1 1 S],1,1); 

% lmiterm([5 1 1 S],1,1); 
% lmiterm([5 1 1 R],-1,1);

LMISYS = getlmis;

% Finding Feasible Solution 
% [~,x] = feasp(LMISYS,[0 1000 1e9 0 0]);
[~,x] = feasp(LMISYS);
R11 = dec2mat(LMISYS,x,1);
R12 = dec2mat(LMISYS,x,2);
R22 = dec2mat(LMISYS,x,3);
R = [R11 R12;R12' R22];
S = dec2mat(LMISYS,x,4);
Gama = dec2mat(LMISYS,x,5);
eig(R-S)

%%

%%
M = eye(np+np)+N'*inv(R)*N;

%%
Q = [R N;N' M];

%%
Ao = blkdiag(Acl,zeros(naw));
Bqo = [Bclq;zeros(naw+np+nc,nu)];
Cyo = [Ccly zeros(nu,naw)];
Dyqo = [Dclyq];
Czo = [Cclz zeros(nz,naw)];
Dzqo = Dclzq;
