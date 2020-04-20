clear 
close all
clc

%% PLANT MODEL
load identified4

load SS_model
A = A_ss_retain_modes_model_1;
B = (B_ss_modal_retain_DOFs_retain_modes(:,3:4))/1e4;
E = B_ss_modal_retain_DOFs_retain_modes(:,1);
C = C_ss_modal_retain_DOFs_retain_modes;

As = A;
Bs = B;
Es = E;
Cs = C;
Gs = G;

% A=A(2:5,2:5);
% B=B(2:5,:);
% E=E(2:5,:);
% C=C(:,2:5);
% G=G(2:5,:);

dum1=size(A);
n=dum1(1,1);
dum2=size(B);
nu=dum2(1,2);
dum3=size(C);
ny=dum3(1,1);
dum4=size(E);
nd=dum4(1,2);

CEp=(((C*E)'*(C*E))^-1)*(C*E)';

H1 = -E*CEp;
H2 = eye(ny,ny)-(C*E)*CEp;
P1 = eye(n,n)-E*CEp*C;
P2 = (eye(ny,ny)-(C*E)*CEp)*C;

%% Finding Y through GA optimization
nvars = n*ny;
lb = -1e1*ones(1,nvars);
ub = 1e1*ones(1,nvars);
PopulationSize_Data = 500;
EliteCount_Data = 10;
CrossoverFraction_Data = 0.9;
MigrationFraction_Data = 0.25;
Generations_Data = 100;
StallGenLimit_Data = 50;
TolFun_Data = 1e-20;
TolCon_Data = 1e-20;
[Y,fval,exitflag,output,population,score] = optimization_for_Y(nvars,lb,ub,PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,MigrationFraction_Data,Generations_Data,StallGenLimit_Data,TolFun_Data,TolCon_Data);

if ny == 2
    Y1 = Y(1,1:6)';
    Y2 = Y(1,7:end)';
%     Y1 = Y(1,1:4)';
%     Y2 = Y(1,5:end)';
    Y = [Y1 Y2];
else
    Y = Y';
end

P = P1+Y*P2;
H = H1+Y*H2;
J = (eye(n,n)+(-E*CEp+Y*(eye(ny,ny)-(C*E)*CEp))*C)*B;

P*E
close all

%% LMI Definition
alfa = 1e1;
beta = 0;

setlmis([]);

X = lmivar(1,[n 1]);
Q2 = lmivar(1,[n 1]);
Th = lmivar(2,[nu n]);
Kh = lmivar(2,[n ny]);
Lh = lmivar(2,[n ny]);
gama_w = lmivar(1,[1 1]);
gama_d = lmivar(1,[1 1]);
gama_v = lmivar(1,[1 1]);

%LMI terms
lmiterm([1 1 1 Th],B,2,'s'); 
lmiterm([1 1 1 X],A,2,'s'); 
lmiterm([1 1 2 0],0);
lmiterm([1 1 3 0],2*G); 
lmiterm([1 1 4 0],2*E);
lmiterm([1 1 5 0],0);
lmiterm([1 1 6 X],1,C');
lmiterm([1 2 2 Q2],1,P*A,'s'); 
lmiterm([1 2 2 Kh],-1,C,'s'); 
lmiterm([1 2 2 0],beta); 
lmiterm([1 2 2 Lh],-1,C,'s'); 
lmiterm([1 2 3 Q2],1,P*G); 
lmiterm([1 2 3 Kh],-1,1); 
lmiterm([1 2 3 Lh],-1,1); 
lmiterm([1 2 4 0],0); 
lmiterm([1 2 5 Q2],1,H); 
lmiterm([1 3 3 gama_w],-1,1);
lmiterm([1 3 4 0],0);
lmiterm([1 3 5 0],0);
lmiterm([1 3 6 0],1);
lmiterm([1 4 4 gama_d],-1,1);
lmiterm([1 5 5 gama_v],-1,1);
lmiterm([1 6 6 0],-alfa);

lmiterm([-2 1 1 X],1,1); 
lmiterm([-2 2 2 X],-1,1); 
lmiterm([-2 2 2 0],0.1); 

lmiterm([-3 1 1 Q2],1,1);

lmiterm([-4 1 1 gama_w],1,1);

lmiterm([-5 1 1 gama_d],1,1);

lmiterm([-6 1 1 gama_v],1,1);

LMISYS = getlmis;

%% Finding Feasible Solution 
n_dec = decnbr(LMISYS);
c = zeros(1,n_dec);
c(1:end) = 1;
c(1:end-1) = 1000;
c(1:end-2) = 3;
options = [1e-6 1000 100 0 0];
% [copt,xopt] = mincx(LMISYS,c,options,[],[]);
% [copt,xopt] = mincx(LMISYS,c);
[copt,xopt] = feasp(LMISYS,options);
% [copt,xopt] = feasp(LMISYS);

X = dec2mat(LMISYS,xopt,1);
Q2 = dec2mat(LMISYS,xopt,2);
Th = dec2mat(LMISYS,xopt,3);
Kh = dec2mat(LMISYS,xopt,4);
Lh = dec2mat(LMISYS,xopt,5);
gama_w = dec2mat(LMISYS,xopt,6)
gama_d = dec2mat(LMISYS,xopt,7)
gama_v = dec2mat(LMISYS,xopt,8)

K = Q2\Kh;
T = Th*X;
N = P*A-K*C;
L = Q2\Lh;

Ak = N+J*T-L*C;
result_opt = [eig(A) eig(N) eig(Ak)]
Bk = L+L*C*H-J*T*H;
Ck = T;
Dk = -T*H;
save('controller','Ak','Bk','Ck','Dk');

Ap = A;
Bp = B;
Cp = C;
Hp = E;
save('plant','Ap','Bp','Cp','Hp');

Hs = Es;
save('plant_s','As','Bs','Cs','Hs');
%%
clear

%%
load plant
load controller
load plant_s

np = sqrt(numel(Ap));
nk = sqrt(numel(Ak));
nu = numel(Bp)/np;
nw = numel(Hp)/np;
ny = numel(Cp)/np;

%%
Plant = ss(Ap,[Bp Hp],Cp,zeros(ny,nu+nw));
K_hat = ss(Ak,[Bk eye(nk) zeros(nk,nu)],Ck,[Dk zeros(nu,nk) eye(nu)]);

systemnames = 'Plant K_hat';
inputvar = '[Nu{2};w{1};Ksi{2+4}]';
outputvar = '[K_hat;Plant;-Nu]';
input_to_Plant = '[K_hat-Nu;w]';
input_to_K_hat = '[Plant;Ksi]';
% cleanupsysic = 'yes';
Ph_SS = sysic;

%%
A = Ph_SS.A;
Bv = Ph_SS.B(:,1:nu);
Bw = Ph_SS.B(:,nu+1:nu+nw);
Bksi = Ph_SS.B(:,nu+nw+1:end);

Cu = Ph_SS.C(1:nu,:);
Duv = Ph_SS.D(1:nu,1:nu);
Duksi = Ph_SS.D(1:nu,nu+nw+1:end);
Duw = Ph_SS.D(1:nu,nu+1:nu+nw);

Cz = Ph_SS.C(nu+1:end,:);
Dzv = Ph_SS.D(nu+1:end,1:nu);
Dzksi = Ph_SS.D(nu+1:end,nu+nw+1:end);
Dzw = Ph_SS.D(nu+1:end,nu+1:nu+nw);

nA = sqrt(numel(A));

%% LMI Definition
% gama = 1e-1;
setlmis([]);

Q = lmivar(1,[nA 1]);
M = lmivar(1,[1 0;1 0]);
delta1 = lmivar(2,[1 1]);
Gama = lmivar(1,[1 1]);
X = lmivar(2,[np+nu nu]);
gama = lmivar(1,[1 0;1 0;1 0;1 0]);

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
lmiterm([1 4 4 gama],-1,1); % LMI #1: -gama*I
lmiterm([1 5 5 delta1],-1,1); % LMI #1: -delta1*I

lmiterm([-2 1 1 Q],1,1); % LMI #3: Q>0

lmiterm([-3 1 1 delta1],1,1); % LMI #4: delta1>0

lmiterm([-4 1 1 Gama],1,1); % LMI #4: Gama>0

lmiterm([-5 1 1 M],1,1); % LMI #4: M>0

lmiterm([-6 1 1 gama],1,1); % LMI #4: M>0

LMISYS = getlmis;

% Finding Feasible Solution 
% [~,x] = feasp(LMISYS,[0 1000 -1 99 0]);
[~,x] = feasp(LMISYS);
Q = dec2mat(LMISYS,x,1);
M = dec2mat(LMISYS,x,2);
delta1 = dec2mat(LMISYS,x,3);
Gama = dec2mat(LMISYS,x,4);
X = dec2mat(LMISYS,x,5);

Lambda = X/M;
K_hat = ss(Ak,[Bk eye(nk) zeros(nk,nu)],Ck,[Dk zeros(nu,nk) eye(nu)]);

Akh = K_hat.A;
Bkh = K_hat.B;
Ckh = K_hat.C;
Dkh = K_hat.D;
save('controller_khat','Akh','Bkh','Ckh','Dkh','Lambda');

eig(Akh)

controller_implemented = ss(Ak,Bk,Ck,Dk);

% set_param('controller_AWBT','AlgebraicLoopSolver','LineSearch')