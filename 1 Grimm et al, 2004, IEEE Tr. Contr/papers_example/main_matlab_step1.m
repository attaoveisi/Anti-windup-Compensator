clear 
close all
clc

%% PLANT MODEL
A = [0 1 0 0;-330.46 -12.15 -2.44 0;0 0 0 1;-812.61 -29.87 -30.1 0];
B = [0 2.71762 0 6.68268]';
H = [0 0 0 15.61]';
C = [1 0 0 0;0 0 1 0];
D = [0 0]';

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = D;
save plant Ap  Bp Hp Cp Dp

%%
nn = 4;

L = [64 2054 -8 -1432;-8 -280 142 10169]';
K = [64.81 213.12 1242.27 85.82];

Ak = A-B*K-L*C;
Bk = L;
Ck = -K;
Dk = zeros(1,2);
save controller.mat Ak Bk Ck Dk nn

%% PLANT MODEL
load controller

Ap = A(1:nn,1:nn);
Bpu = B(1:nn,:);
Cpy = C(:,1:nn);
Cpz = [0 0 1 0];
Bpw = H(1:nn,:);

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
controller = ss(Ak,[zeros(np,nw),Bk,[eye(nk) zeros(nk,nu)]],Ck,[zeros(nu,nw),Dk,[zeros(nu,nk) eye(nu)]]);

systemnames = 'plant controller';
inputvar = '[w;v(5)]';
outputvar = '[plant(1);controller]';
input_to_plant = '[w;controller]';
input_to_controller = '[w;plant(2:3);v]';
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