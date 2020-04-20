clear 
close all
clc

%% PLANT MODEL
load identified4

As = A;
Bs = B;
Es = E;
Cs = C;
Gs = G;

A=A(2:5,2:5);
B=B(2:5,:);
E=E(2:5,:);
C=C(:,2:5);
G=G(2:5,:);

Ap = A;
Bp = B;
Cp = C;
Hp = E;

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
lb = -1e2*ones(1,nvars);
ub = 1e2*ones(1,nvars);
PopulationSize_Data = 500;
EliteCount_Data = 10;
CrossoverFraction_Data = 0.9;
MigrationFraction_Data = 0.25;
Generations_Data = 200;
StallGenLimit_Data = 60;
TolFun_Data = 1e-30;
TolCon_Data = 1e-30;
[Y,fval,exitflag,output,population,score] = optimization_for_Y1(nvars,lb,ub,PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,MigrationFraction_Data,Generations_Data,StallGenLimit_Data,TolFun_Data,TolCon_Data);

if ny == 2
%     Y1 = Y(1,1:6)';
%     Y2 = Y(1,7:end)';
    Y1 = Y(1,1:4)';
    Y2 = Y(1,5:end)';
    Y = [Y1 Y2];
else
    Y = Y';
end
% load Y
P = P1+Y*P2;
H = H1+Y*H2;
J = (eye(n,n)+(-E*CEp+Y*(eye(ny,ny)-(C*E)*CEp))*C)*B;

P*E
close all

%% LMI Definition
setlmis([]);

Q2 = lmivar(1,[n 1]);
Kh = lmivar(2,[n ny]);
Lh = lmivar(2,[n ny]);
gama_w = lmivar(1,[1 1]);
gama_v = lmivar(1,[1 1]);

%LMI terms
lmiterm([1 1 1 Q2],1,P*A,'s'); 
lmiterm([1 1 1 Kh],-1,C,'s'); 
lmiterm([1 1 1 Lh],-1,C,'s'); 
lmiterm([1 1 1 0],1); 
lmiterm([1 1 2 Q2],1,P*G); 
lmiterm([1 1 2 Kh],-1,1); 
lmiterm([1 1 2 Lh],-1,1); 
lmiterm([1 1 3 Q2],1,H); 
lmiterm([1 2 2 gama_w],-1,1);
lmiterm([1 3 3 gama_v],-1,1);
lmiterm([-2 1 1 Q2],1,1);
lmiterm([-3 1 1 gama_w],1,1);
lmiterm([-4 1 1 gama_v],1,1);

LMISYS = getlmis;

options = [1e-12 1000 50 0 0];
% [copt,xopt] = feasp(LMISYS,options);
[copt,xopt] = feasp(LMISYS);

Q2 = dec2mat(LMISYS,xopt,1);
Kh = dec2mat(LMISYS,xopt,2);
Lh = dec2mat(LMISYS,xopt,3);
gama_w = dec2mat(LMISYS,xopt,4);
gama_v = dec2mat(LMISYS,xopt,5);

K = Q2\Kh;
N = P*A-K*C;
L = Q2\Lh;

save import2 n nu ny nd Y P H J A B E C G K N L Q2 Kh Lh