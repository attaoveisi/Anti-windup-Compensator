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

%% PLANT MODEL
Dpw = zeros(ny,nw);
plant = ss(Ap,[Hp,Bp],Cp,[Dpw,Dp]);
controller = ss(Ak,[zeros(nk,nw),Bk],Ck,[zeros(nu,nw),Dk]);

systemnames = 'plant controller';
inputvar = '[w]';
outputvar = '[controller]';
input_to_plant = '[w;controller]';
input_to_controller = '[w;plant]';
cleanupsysic = 'yes';
sys_cl = sysic;
                                                                                                                                                                      
Ad = sys_cl.A;
nd = sqrt(numel(Ad));
Bd = sys_cl.B;
Cd = sys_cl.C;
Dd = sys_cl.D;

% Ky = tf(rlqg);
% Kw = [0;0];
% Pu = tf(ss(Ap,Bp,Cp,Dp));
% Pw = tf(ss(Ap,Hp,Cp,Dpw));
% Pd = ((eye(2)-Ky*Pu)^(-1))*(Ky*Pw);
% figure;
% bode(Pd)
% hold on
% bode(ss(Ad,Bd,Cd,Dd))

%%
Ao = blkdiag(Ap,Ad);
WA = [zeros(np,np),Bp*Cd;zeros(nd,np),zeros(nd,nd)];
WB = [Bp*Dd,zeros(np,ny);Bd,zeros(nd,ny)];
WC = [zeros(np,nw),Cp';zeros(nd,nw),Cd'*Dp'];

%% LMI
% Define variables
P11 = sdpvar(np,np);
P12 = sdpvar(np,nd);
P22 = sdpvar(nd,nd);
gamad = sdpvar(1);

P = [P11,P12;P12',P22];
WD = [-gamad*eye(nw,nw),Dd'*Dp';Dp*Dd,-gamad*eye(ny,ny)];

% Define constraints and objective
Pi1 = [P*Ao+Ao'*P+P*WA+WA'*P,WC+P*WB;WC'+WB'*P,WD];
Pi2 = [Ad'*P22+P22*Ad,P22*Bd;Bd'*P22,-gamad*eye(nw,nw)];

% ,P(:) >= -30,P(:) <= 30
Constraints = [Pi1 <= 0,Pi2 <= 0,gamad >= 0,P >= 0.001*eye(size(P)),P22(:) >= -1,P22(:) <= 1,P11(:) >= -1,P11(:) <= 1];
% Constraints = [Pi1 <= 0,Pi2 <= 0,gamad >= 0,P >= 0];
Objective = [];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');
% options = sdpsettings('verbose',1,'solver','mosek','radius',2.5);

optimize(Constraints,Objective,options);
P11 = value(P11);
P12 = value(P12);
P22 = value(P22);
gamad = value(gamad)

P = [P11,P12;P12',P22];
WD = [-gamad*eye(nw,nw),Dd'*Dp';Dp*Dd,-gamad*eye(ny,ny)];
Pi1 = [P*Ao+Ao'*P+P*WA+WA'*P,WC+P*WB;WC'+WB'*P,WD];
Pi2 = [Ad'*P22+P22*Ad,P22*Bd;Bd'*P22,-gamad*eye(nw,nw)];

eig(P)
eig(Pi1)
eig(Pi2)

%%
Bo = [zeros(np,nu),zeros(np,nw);zeros(nd,nu),Bd];
Cdo = [zeros(nu,np),Cd;zeros(nw,np),zeros(nw,nd)];
Ddo = [zeros(nu,nu),Dd;zeros(nw,nu),zeros(nw,nw)];
Cpo = [Cp zeros(ny,nd)];
% G1 = [eye(np,np),zeros(np,nd);zeros(nu,np),zeros(nu,nd)];
W = eye(nu,nu);
Wt = blkdiag(W,eye(nw,nw));
Itnw = eye(nu+nw,nu+nw);

Sai = [Ao'*P+P*Ao,P*Bo+Cdo'*Wt,Cpo';Bo'*P+Wt*Cdo,Wt*Ddo+Ddo'*Wt-gamad*Itnw,zeros(nu+nw,ny);Cpo,zeros(ny,nu+nw),-gamad*eye(ny,ny)];
H = [Bp' zeros(nu,nd) -eye(nu,nu) zeros(nu,nw) Dp']*blkdiag(P,Wt,eye(ny,ny));
G = [eye(np,np) zeros(np,nd) zeros(np,nu) zeros(np,nw) zeros(np,ny);zeros(nu,np) zeros(nu,nd) eye(nu,nu) zeros(nu,nw) zeros(nu,ny)];

%%
% Define variables
% E_dum = sdpvar(nu,1);
% E = diag(E_dum);
E = sdpvar(nu,nu,'full');
F = sdpvar(nu,np);
GAMA = [F E];

% Define constraints and objective
Pi1 = Sai+H'*GAMA*G+G'*GAMA'*H;

aa = 100;
% , GAMA(:) >= -aa,  aa >= GAMA(:)
% ,E(:) >= -aa,  aa >= E(:)

Constraints = [Pi1 <= 0];
Objective = [];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');
% options = sdpsettings('verbose',1,'solver','mosek','radius',4000);

optimize(Constraints,Objective,options);
% E_dum = value(E_dum);
% E = diag(E_dum);
E = value(E)*0.0001;
F = value(F);
GAMA = [F E];

% AW compensator structure
dumsys1 = ss(Ap+Bp*F,Bp*E,F,E-eye(size(E)));
dumsys2 = ss(Ap+Bp*F,Bp*E,Cp+Dp*F,Dp*E);
Ms = (simplify(tf(dumsys1)-eye(size(E))))*inv(E);
Ns = tf(dumsys2)*inv(E);

Pi1 = Sai+H'*GAMA*G+G'*GAMA'*H;
eig(Pi1)

[eig(Ap) eig(Ap+Bp*F)]

%% LMI Definition
setlmis([]);

GAMA = lmivar(2,[nu nu+np]);

%LMI terms
lmiterm([1 1 1 GAMA],H',G,'s'); 
lmiterm([1 1 1 0],Sai);

LMISYS = getlmis;

[~,x] = feasp(LMISYS,[0 1000 500 99 0]);
% [~,x] = feasp(LMISYS);
GAMA = dec2mat(LMISYS,x,1);
F = GAMA(:,1:np);
E = GAMA(:,np+1:end)*0.1;
% E = eye(nu);

dumsys1 = ss(Ap+Bp*F,Bp*E,F,E-eye(size(E)));
dumsys2 = ss(Ap+Bp*F,Bp*E,Cp+Dp*F,Dp*E);

Pi1 = Sai+H'*GAMA*G+G'*GAMA'*H;
eig(Pi1)

[eig(Ap) eig(Ap+Bp*F)]

    F1 = tf(400*2*pi,[1 400*2*pi]);
    F2 = tf(100*2*pi,[1 100*2*pi]);