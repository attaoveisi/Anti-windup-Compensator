clc
clear
close all

%%
load plant
load controller

np = sqrt(numel(Ap));
nk = sqrt(numel(Ak));
nu = numel(Bp)/np;

%%
load Mathematica2MatlabNum
load Mathematica2MatlabDenum

Ph_TF_num = {PhTFnum11,PhTFnum12,PhTFnum13,PhTFnum14,PhTFnum15,PhTFnum16,PhTFnum17,PhTFnum18,PhTFnum19,PhTFnum110,PhTFnum111;
             PhTFnum21,PhTFnum22,PhTFnum23,PhTFnum24,PhTFnum25,PhTFnum26,PhTFnum27,PhTFnum28,PhTFnum29,PhTFnum210,PhTFnum211;
             PhTFnum31,PhTFnum32,PhTFnum33,PhTFnum34,PhTFnum35,PhTFnum36,PhTFnum37,PhTFnum38,PhTFnum39,PhTFnum310,PhTFnum311;
            PhTFnum41,PhTFnum42,PhTFnum43,PhTFnum44,PhTFnum45,PhTFnum46,PhTFnum47,PhTFnum48,PhTFnum49,PhTFnum410,PhTFnum411;
            PhTFnum51,PhTFnum52,PhTFnum53,PhTFnum54,PhTFnum55,PhTFnum56,PhTFnum57,PhTFnum58,PhTFnum59,PhTFnum510,PhTFnum511;
            PhTFnum61,PhTFnum62,PhTFnum63,PhTFnum64,PhTFnum65,PhTFnum66,PhTFnum67,PhTFnum68,PhTFnum69,PhTFnum610,PhTFnum611};
        
Ph_TF_denum = {PhTFdenum11,PhTFdenum12,PhTFdenum13,PhTFdenum14,PhTFdenum15,PhTFdenum16,PhTFdenum17,PhTFdenum18,PhTFdenum19,PhTFdenum110,PhTFdenum111;
               PhTFdenum21,PhTFdenum22,PhTFdenum23,PhTFdenum24,PhTFdenum25,PhTFdenum26,PhTFdenum27,PhTFdenum28,PhTFdenum29,PhTFdenum210,PhTFdenum211;
               PhTFdenum31,PhTFdenum32,PhTFdenum33,PhTFdenum34,PhTFdenum35,PhTFdenum36,PhTFdenum37,PhTFdenum38,PhTFdenum39,PhTFdenum310,PhTFdenum311;
               PhTFdenum41,PhTFdenum42,PhTFdenum43,PhTFdenum44,PhTFdenum45,PhTFdenum46,PhTFdenum47,PhTFdenum48,PhTFdenum49,PhTFdenum410,PhTFdenum411;
               PhTFdenum51,PhTFdenum52,PhTFdenum53,PhTFdenum54,PhTFdenum55,PhTFdenum56,PhTFdenum57,PhTFdenum58,PhTFdenum59,PhTFdenum510,PhTFdenum511;
               PhTFdenum61,PhTFdenum62,PhTFdenum63,PhTFdenum64,PhTFdenum65,PhTFdenum66,PhTFdenum67,PhTFdenum68,PhTFdenum69,PhTFdenum610,PhTFdenum611};

%%
Ph_TF = tf(Ph_TF_num,Ph_TF_denum);
Ph_SS = ss(Ph_TF,'minimal');
% hsvd(Ph_SS)

%%
A = Ph_SS.A;
Bv = Ph_SS.B(:,1:2);
Bw = Ph_SS.B(:,3);
Bksi = Ph_SS.B(:,4:end);

Cu = Ph_SS.C(1:2,:);
Duv = Ph_SS.D(1:2,1:2);
Duksi = Ph_SS.D(1:2,4:end);
Duw = Ph_SS.D(1:2,3);

Cz = Ph_SS.C(3,:);
Dzv = Ph_SS.D(3,1:2);
Dzksi = Ph_SS.D(3,4:end);
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