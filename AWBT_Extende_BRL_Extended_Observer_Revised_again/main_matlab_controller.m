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

lmiterm([-6 1 1 gama],1,1); % LMI #4: gama>0

lmiterm([7 1 1 gama],1,1);
lmiterm([7 1 1 0],-2e6);

LMISYS = getlmis;

% Finding Feasible Solution 
% [~,x] = feasp(LMISYS,[0 1000 1e9 0 0]);
[~,x] = feasp(LMISYS);
Q = dec2mat(LMISYS,x,1);
M = dec2mat(LMISYS,x,2);
delta1 = dec2mat(LMISYS,x,3);
Gama = dec2mat(LMISYS,x,4);
X = dec2mat(LMISYS,x,5);
gama = dec2mat(LMISYS,x,6);

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