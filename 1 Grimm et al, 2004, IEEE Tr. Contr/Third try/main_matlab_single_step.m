clear 
close all
clc

%% PLANT MODEL
% load identified
% load identified_freq
load identified_t1
load identified_t1_full
% B = B/1e4;

% A = A(1:4,1:4);
% B = B(1:4,:);
% H = H(1:4,:);
% C = C(:,1:4);

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = zeros(1,2);
save plant Ap Bp Hp Cp Dp

%%
Bw=H;
sys1=ss(A,Bw,C,0);

%%%%%%% Kalma Filter Design
sysk=ss(A,[B Bw],C,0);

[kest,L,P] = kalman(sysk,5e3,1,0);

%%%%%%% Optimal Gain
syst=ss(A,B,C,0);
[K,~,e] = lqry(syst,1e4,1e0*eye(2),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ak = rlqg.A;
Bk = rlqg.B;
Ck = rlqg.C;
Dk = rlqg.D;
save controller.mat Ak Bk Ck Dk

%% PLANT MODEL
% load plant
% load controller

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

plant = ss(Ap,[Bpw,Bpu],[Cpz;Cpy],[Dpzw,Dpzu;Dpyw,Dpyu]);
controller = ss(Ak,[zeros(nk,nw),Bk,[eye(nk) zeros(nk,nu)]],Ck,[zeros(nu,nw),Dk,[zeros(nu,nk) eye(nu)]]);

systemnames = 'plant controller';
inputvar = '[w{1};v{14}]';
outputvar = '[plant(1);controller(1:2)]';
input_to_plant = '[w(1);controller(1:2)]';
input_to_controller = '[w(1);plant(2);v(1:14)]';
cleanupsysic = 'yes';
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

%%
clc
clear
load('plant.mat');
load('controller.mat');
load('closed_loop.mat');

Bpu = Bp;
Cpy = Cp;
Cpz = Cpy;
Bpw = Hp;

dum1=size(Ap);
np=dum1(1,1);
dum2=size(Bpu);
nu=dum2(1,2);
dum3=size(Cpy);
ny=dum3(1,1);
nz=ny;
dum4=size(Bpw);
nw=dum4(1,2);

dum5=size(Ak);
nc=dum5(1,2);

naw=np+nc;
ncl=np+nc;

Dpyu = zeros(ny,nu);
Dpyw = zeros(ny,nw);
Dpzu = zeros(nz,nu);
Dpzw = zeros(nz,nw);

% Define variables
R11 = sdpvar(np,np);
R12 = sdpvar(np,nc);
R22 = sdpvar(nc,nc);
S = sdpvar(np+nc,np+nc);
gama = sdpvar(1);

% Define constraints and objective
Pi1 = [R11*Ap'+Ap*R11,Bpw,R11*Cpz';Bpw',-gama*eye(nw,nw),Dpzw';Cpz*R11,Dpzw,-gama*eye(nz,nz)];
Pi2 = [S*Acl'+Acl*S,Bclw,S*Cclz';Bclw',-gama*eye(nw,nw),Dclzw';Cclz*S,Dclzw,-gama*eye(nz,nz)];
Pi3 = [R11,R12;R12',R22];
Pi4 = ([R11,R12;R12',R22])-S;

Constraints = [Pi1 <= -.1*eye(size(Pi1)),Pi2 <= 0,Pi3 >= 0,Pi4 >= 0,gama >= 0];
% Objective = [];
Objective = gama;

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

% Solve the problem

%%%%%%%%%% Method 1 %%%%%%%%%%%

% sol = optimize(Constraints,Objective,options);
% Analyze error flags
% if sol.problem == 0
%  R11 = value(R11);
%  R12 = value(R12);
%  R22 = value(R22);
%  S = value(S);
%  gama = value(gama);
% else
%  sol.info
% end

%%%%%%%%%% Method 2 %%%%%%%%%%%

optimize(Constraints,gama,options);
R11 = value(R11);
R12 = value(R12);
R22 = value(R22);
S = value(S);
gama = value(gama);

%% Test the results
Pi1 = [R11*Ap'+Ap*R11,Bpw,R11*Cpz';Bpw',-gama*eye(nw,nw),Dpzw';Cpz*R11,Dpzw,-gama*eye(nz,nz)];
Pi2 = [S*Acl'+Acl*S,Bclw,S*Cclz';Bclw',-gama*eye(nw,nw),Dclzw';Cclz*S,Dclzw,-gama*eye(nz,nz)];
Pi3 = [R11,R12;R12',R22];
Pi4 = ([R11,R12;R12',R22])-S;

eig(Pi1)
eig(Pi2)
eig(Pi3)
eig(Pi4)
gama

%%
R = [R11,R12;R12',R22];
if naw == np+nc  
    N = chol(R*inv(S)*R-R);
    N = N';
else
    [N,Ndum] = fullrf(R*inv(S)*R-R);
end

max(max(abs(R*inv(S)*R-R-N*N')))/min(min(abs(R*inv(S)*R-R)))

M = eye(naw,naw)+N'*R^-1*N;
Q = [R,N;N',M];

n = naw+nc+np;
nv = nc+nu;

Ao = blkdiag(Acl,zeros(naw,naw));
Bqo = [Bclq;zeros(naw,nu)];
Cyo = [Ccly,zeros(nu,naw)];
Dyqo = Dclyq;
Czo = [Cclz,zeros(nz,naw)];
Dzqo = Dclzq;
H1 = [zeros(np+nc,naw),Bclv;eye(naw,naw),zeros(naw,nc+nu)]';
G1 = zeros(naw+nu,n);
G1(1:naw,end-naw+1:end) = eye(naw,naw);
G2 = [zeros(naw,nu);eye(nu,nu)];
H2 = [zeros(nu,naw),Dclyv]';
H3 = [zeros(nz,naw),Dclzv]';
Bw = [Bclw;zeros(naw,nw)];
Dzw = Dclzw;
Dyw = Dclyw;

delta = 1e-1;
Vs = 1e0*eye(nu,nu);
U = delta*inv(Vs);
Sai = [Q*Ao'+Ao*Q,Bqo*U+Q*Cyo',Bw,Q*Czo';U*Bqo'+Cyo*Q,Dyqo*U+U*Dyqo'-2*U,Dyw,U*Dzqo';Bw',Dyw',-gama*eye(nw,nw),Dzw';Czo*Q,Dzqo*U,Dzw,-gama*eye(nz,nz)];
H = [H1,H2,zeros(naw+nv,nw),H3];
G = [G1*Q,G2*U,zeros(naw+nu,nw+nz)];

%% Step 3
nu = 2;

delta = 1e3;
Vs = 1e0*eye(nu,nu);
U = delta*inv(Vs);
Sai = [Q*Ao'+Ao*Q,Bqo*U+Q*Cyo',Bw,Q*Czo';U*Bqo'+Cyo*Q,Dyqo*U+U*Dyqo'-2*U,Dyw,U*Dzqo';Bw',Dyw',-gama*eye(nw,nw),Dzw';Czo*Q,Dzqo*U,Dzw,-gama*eye(nz,nz)];
H = [H1,H2,zeros(naw+nv,nw),H3];
G = [G1*Q,G2*U,zeros(naw+nu,nw+nz)];

%% LMI Definition
setlmis([]);

GAMA = lmivar(2,[naw+nv naw+nu]);

%LMI terms
lmiterm([1 1 1 GAMA],H',G,'s'); % LMI #1:
lmiterm([1 1 1 0],Sai); % LMI #1:

LMISYS = getlmis;

%% Define variables
GAMA1 = sdpvar(naw,naw,'full');
GAMA2 = sdpvar(naw,nu);
GAMA3 = sdpvar(nv,naw);
GAMA4 = sdpvar(nv,nu);
% GAMA4 = zeros(nv,nu);

GAMA = [GAMA1,GAMA2;GAMA3,GAMA4];

% Define constraints and objective
Constraints = [H'*GAMA*G+G'*GAMA'*H+Sai <= 0];

% Solve the problem
optimize(Constraints);
GAMA1 = value(GAMA1);
GAMA2 = value(GAMA2);
GAMA3 = value(GAMA3);
GAMA4 = value(GAMA4);
GAMA = [GAMA1,GAMA2;GAMA3,GAMA4];

eig(H'*GAMA*G+G'*GAMA'*H+Sai)
load identified_t1_full
sim('AW_test')

%% Finding Feasible Solution 
% [~,x] = feasp(LMISYS);
% [~,x] = feasp(LMISYS,[0 1000 9e3 0 0]);
% GAMA = dec2mat(LMISYS,x,1);
% GAMA1 = GAMA(1:naw,1:naw);
% GAMA2 = GAMA(1:naw,naw+1:end);
% GAMA3 = GAMA(naw+1:end,1:naw);
% GAMA4 = GAMA(naw+1:end,naw+1:end);
% eig(GAMA1)
% GAMA0 = GAMA;
% save GAMA GAMA1 GAMA2 GAMA3 GAMA4
% GAMA10 = GAMA1;
% GAMA20 = GAMA2;
% GAMA30 = GAMA3;
% GAMA40 = GAMA4;
% save GAMA0 GAMA10 GAMA20 GAMA30 GAMA40
% dumsys = ss(GAMA1,GAMA2,GAMA3,GAMA4);
% 
load identified_t1_full