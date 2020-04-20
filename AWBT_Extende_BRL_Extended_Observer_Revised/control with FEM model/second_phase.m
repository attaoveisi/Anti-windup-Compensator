clc
clear
close all

%%
load plant
load controller
load additional

np = sqrt(numel(Ap));
nk = sqrt(numel(Ak));
nu = numel(Bp)/np;

%%
syms s
P11h = -((eye(2)-((Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp)))^(-1))*(Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp);
P12h = +((eye(2)-((Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp)))^(-1))*(Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Hp);
P13h = +((eye(2)-((Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp)))^(-1))*(Ck*(s*eye(nk)-Ak)^(-1));
P14h = +((eye(2)-((Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp)))^(-1));
P21h = (Cp*(s*eye(np)-Ap)^(-1)*Bp)*P11h-(Cp*(s*eye(np)-Ap)^(-1)*Bp);
P22h = (Cp*(s*eye(np)-Ap)^(-1)*Bp)*P12h-(Cp*(s*eye(np)-Ap)^(-1)*Hp);
P23h = (Cp*(s*eye(np)-Ap)^(-1)*Bp)*P13h;
P24h = (Cp*(s*eye(np)-Ap)^(-1)*Bp)*P14h;
P31h = -eye(2);
P32h = zeros(2,1);
P33h = zeros(2,4);
P34h = zeros(2,2);

Ph = [P11h P12h P13h P14h;
      P21h P22h P23h P24h;
      P31h P32h P33h P34h];
[nPh,mPh] = size(Ph);

for i = 1:nPh
    for j = 1:mPh
        Ph_TF_num{i,j} = syms2num(Ph(i,j));
        Ph_TF_denum{i,j} = syms2denum(Ph(i,j));
    end
end

Ph_TF = tf(Ph_TF_num,Ph_TF_denum);
Ph_SS = ss(Ph_TF,'minimal');

A = Ph_SS.A;
Bv = Ph_SS.B(:,1:2);
Bw = Ph_SS.B(:,3);
Bksi = Ph_SS.B(:,4:end);
% Bksi = -Ph_SS.B(:,4:end);

Cu = Ph_SS.C(1:2,:);
Duv = Ph_SS.D(1:2,1:2);
Duksi = Ph_SS.D(1:2,4:end);
% Duksi = -Ph_SS.D(1:2,4:end);
Duw = Ph_SS.D(1:2,3);

Cz = Ph_SS.C(3,:);
Dzv = Ph_SS.D(3,1:2);
Dzksi = Ph_SS.D(3,4:end);
% Dzksi = -Ph_SS.D(3,4:end);
Dzw = Ph_SS.D(3,3);

nA = sqrt(numel(A));

%% LMI Definition
gama = 1e1;

setlmis([]);

Q = lmivar(1,[nA 1]);
M = lmivar(1,[1 0;1 0]);
delta1 = lmivar(2,[1 1]);
Gama = lmivar(1,[1 1]);
X = lmivar(2,[np+nu nu]);

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

lmiterm([-5 1 1 M],1,1); % LMI #4: M>0

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
Khat = ss(Ak,[Bk eye(4) zeros(4,2)],Ck,[Dk zeros(2,4) eye(2)]);
% Khat = ss(Ak,[Bk eye(4)],Ck,[Dk zeros(2,4)]);

Akh = Khat.A;
Bkh = Khat.B;
Ckh = Khat.C;
Dkh = Khat.D;
save('controller_khat','Akh','Bkh','Ckh','Dkh','Lambda');

controller_implemented = ss(Ak,Bk,Ck,Dk);

% set_param('controller_AWBT','AlgebraicLoopSolver','LineSearch')

%% Kalman Filter
%%%%%%% Kalma Filter Design
sysk=ss(Ap,[Bp Hp],Cp,0);
[kest,Lest,Pest] = kalman(sysk,1e4,1,0);

%%%%%%% Optimal Gain
syst=ss(Ap,Bp,Cp,0);
[KLQR,SLQR,eLQR] = lqry(syst,4e2,1e0*eye(2),0);

%%%%%%% Controller
rlqg = lqgreg(kest,KLQR);

AkLQG = rlqg.A;
BkLQG = rlqg.B;
CkLQG = rlqg.C;
DkLQG = rlqg.D;
save controller_LQG.mat AkLQG BkLQG CkLQG DkLQG

%% PLANT MODEL
load plant
load controller_LQG

%%
syms s
Ak = AkLQG;
Bk = BkLQG;
Ck = CkLQG;
Dk = DkLQG;
P11h = -((eye(2)-((Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp)))^(-1))*(Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp);
P12h = +((eye(2)-((Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp)))^(-1))*(Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Hp);
P13h = +((eye(2)-((Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp)))^(-1))*(Ck*(s*eye(nk)-Ak)^(-1));
P14h = +((eye(2)-((Ck*(s*eye(nk)-Ak)^(-1)*Bk+Dk)*(Cp*(s*eye(np)-Ap)^(-1)*Bp)))^(-1));
P21h = (Cp*(s*eye(np)-Ap)^(-1)*Bp)*P11h-(Cp*(s*eye(np)-Ap)^(-1)*Bp);
P22h = (Cp*(s*eye(np)-Ap)^(-1)*Bp)*P12h-(Cp*(s*eye(np)-Ap)^(-1)*Hp);
P23h = (Cp*(s*eye(np)-Ap)^(-1)*Bp)*P13h;
P24h = (Cp*(s*eye(np)-Ap)^(-1)*Bp)*P14h;
P31h = -eye(2);
P32h = zeros(2,1);
P33h = zeros(2,4);
P34h = zeros(2,2);

Ph = [P11h P12h P13h P14h;
      P21h P22h P23h P24h;
      P31h P32h P33h P34h];
[nPh,mPh] = size(Ph);

for i = 1:nPh
    for j = 1:mPh
        Ph_TF_num{i,j} = syms2num(Ph(i,j));
        Ph_TF_denum{i,j} = syms2denum(Ph(i,j));
    end
end

Ph_TF = tf(Ph_TF_num,Ph_TF_denum);
Ph_SS = ss(Ph_TF,'minimal');

A = Ph_SS.A;
Bv = Ph_SS.B(:,1:2);
Bw = Ph_SS.B(:,3);
Bksi = Ph_SS.B(:,4:end);
% Bksi = -Ph_SS.B(:,4:end);

Cu = Ph_SS.C(1:2,:);
Duv = Ph_SS.D(1:2,1:2);
Duksi = Ph_SS.D(1:2,4:end);
% Duksi = -Ph_SS.D(1:2,4:end);
Duw = Ph_SS.D(1:2,3);

Cz = Ph_SS.C(3,:);
Dzv = Ph_SS.D(3,1:2);
Dzksi = Ph_SS.D(3,4:end);
% Dzksi = -Ph_SS.D(3,4:end);
Dzw = Ph_SS.D(3,3);

nA = sqrt(numel(A));

%% LMI Definition
gama = 1e-2;

setlmis([]);

Q = lmivar(1,[nA 1]);
M = lmivar(1,[1 0;1 0]);
delta1 = lmivar(2,[1 1]);
Gama = lmivar(1,[1 1]);
X = lmivar(2,[np+nu nu]);

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

lmiterm([-5 1 1 M],1,1); % LMI #4: M>0

LMISYS = getlmis;

%% Finding Feasible Solution 
% [~,x] = feasp(LMISYS,[0 1000 -1 99 0]);
[~,x] = feasp(LMISYS);
Q = dec2mat(LMISYS,x,1);
M = dec2mat(LMISYS,x,2);
delta1 = dec2mat(LMISYS,x,3);
Gama = dec2mat(LMISYS,x,4);
X = dec2mat(LMISYS,x,5);

Lambda_LQG = X/M;
Khat_LQG = ss(Ak,[Bk eye(4) zeros(4,2)],Ck,[Dk zeros(2,4) eye(2)]);
