clear 
close all
clc

%% PLANT MODEL
load identified

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = D;
save plant Ap  Bp Hp Cp Dp

%%
nn = 16;

% A = A(1:nn,1:nn);
% B = B(1:nn,:);
% C = C(:,1:nn);
% H = H(1:nn,:);

nnn = nn;

A = A(nn-nnn+1:nn,nn-nnn+1:nn);
B = B(nn-nnn+1:nn,:);
C = C(:,nn-nnn+1:nn);
H = H(nn-nnn+1:nn,:);

Bw=H;
sys1=ss(A,Bw,C,0);

%%%%%%% Kalma Filter Design
sysk=ss(A,[B Bw],C,0);

[kest,L,P] = kalman(sysk,1e3,1,0);

%%%%%%% Optimal Gain
syst=ss(A,B,C,0);
[K,S,e] = lqry(syst,3.5e0,1e0*eye(2),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ak = rlqg.A;
Bk = rlqg.B;
Ck = rlqg.C;
Dk = rlqg.D;
save controller.mat Ak Bk Ck Dk nn

%% PLANT MODEL
load plant
load controller

% Ap = Ap(1:nn,1:nn);
% Bpu = Bp(1:nn,:);
% Cpy = Cp(:,1:nn);
% Cpz = Cpy;
% Bpw = Hp(1:nn,:);

Ap = Ap(nn-nnn+1:nn,nn-nnn+1:nn);
Bpu = Bp(nn-nnn+1:nn,:);
Cpy = Cp(:,nn-nnn+1:nn);
Cpz = Cpy;
Bpw = Hp(nn-nnn+1:nn,:);


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

plant = ss(Ap,[Bpw,Bpu],[Cpz;Cpy],[Dpzw,Dpzu;Dpyw,Dpyu]);
controller = ss(Ak,[zeros(nk,nw),Bk,[eye(nk) zeros(nk,nu)]],Ck,[zeros(nu,nw),Dk,[zeros(nu,nk) eye(nu)]]);

systemnames = 'plant controller';
inputvar = '[w;v(18)]';
outputvar = '[plant(1);controller]';
input_to_plant = '[w;controller]';
input_to_controller = '[w;plant(2);v]';
% cleanupsysic = 'yes';
sys_cl = sysic;
                                                                                                                                                                      
Acl = sys_cl.A;
ncl = sqrt(numel(Acl));
Bclw = sys_cl.B(:,1:nw);
Bclv = sys_cl.B(:,nw+1:end);
Bclq = zeros(ncl,nu);
Cclz = sys_cl.C(1:nz,:);
Ccly = sys_cl.C(nz+1:end,:);
Dclzw = sys_cl.D(1:nz,1:nw);
Dclzv = sys_cl.D(1:nz,nw+1:end);
Dclyw = sys_cl.D(nz+1:end,1:nw);
Dclyv = sys_cl.D(nz+1:end,nw+1:end);
Dclzq = zeros(nz,nu);
Dclyq = zeros(nu,nu);


save closed_loop.mat Acl Bclw Bclv Bclq Cclz Ccly Dclzw Dclzv Dclzq Dclyw Dclyv Dclyq