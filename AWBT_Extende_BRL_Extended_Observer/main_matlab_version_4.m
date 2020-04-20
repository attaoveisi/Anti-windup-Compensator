clear 
close all
clc

%% PLANT MODEL
load identified2

A=A(1:4,1:4);
B=B(1:4,:);
E=E(1:4,:);
G=G(1:4,:);
C=C(:,1:4);

dum1=size(A);
n=dum1(1,1);
dum2=size(B);
nu=dum2(1,2);
dum3=size(C);
ny=dum3(1,1);
dum4=size(E);
nd=dum4(1,2);

CEp=((C*E)'*(C*E))^-1;

H1 = -E*CEp;
H2 = eye(ny,ny)-(C*E)*CEp;
P1 = eye(n,n)-E*CEp*C;
P2 = (eye(ny,ny)-(C*E)*CEp)*C;

%% Finding Y through GA optimization
nvars = n;
lb = -1e1*ones(1,n);
ub = 1e1*ones(1,n);
PopulationSize_Data = 500;
EliteCount_Data = 10;
CrossoverFraction_Data = 0.9;
MigrationFraction_Data = 0.25;
Generations_Data = 2000;
StallGenLimit_Data = 100;
TolFun_Data = 1e-20;
TolCon_Data = 1e-20;
[Y,fval,exitflag,output,population,score] = optimization_for_Y(nvars,lb,ub,PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,MigrationFraction_Data,Generations_Data,StallGenLimit_Data,TolFun_Data,TolCon_Data);

Y = Y';
P = P1+Y*P2;
H = H1+Y*H2;
J = (eye(n,n)+(-E*CEp+Y*(eye(ny,ny)-(C*E)*CEp))*C)*B;

P*E
close all

%% LMI Definition
alfa = 1;
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

lmiterm([-3 1 1 Q2],1,1);

lmiterm([-4 1 1 gama_w],1,1);

lmiterm([-5 1 1 gama_d],1,1);

lmiterm([-6 1 1 gama_v],1,1);

LMISYS = getlmis;

%% Finding Feasible Solution 
n_dec = decnbr(LMISYS);
c = zeros(1,n_dec);
c(1:end) = -10;
c(1:end-1) = -10;
c(1:end-2) = -1e-3;
options = [1e-6 1000 0 0 0];
[copt,xopt] = mincx(LMISYS,c,options,[],[]);

X = dec2mat(LMISYS,xopt,1);
Q2 = dec2mat(LMISYS,xopt,2);
Th = dec2mat(LMISYS,xopt,3);
Kh = dec2mat(LMISYS,xopt,4);
Lh = dec2mat(LMISYS,xopt,5);
gama_w = dec2mat(LMISYS,xopt,6)
gama_d = dec2mat(LMISYS,xopt,7)
gama_v = dec2mat(LMISYS,xopt,8)

K = Q2\Kh;
T = Th/X;
N = P*A-K*C;
L = Q2\Lh;

result_opt = [eig(A) eig(A+B*T) eig(N)]

Ak = N-J*T-L*C;
eig(Ak)
Bk = L+L*C*H-J*T*H;
Ck = T;
Dk = -T*H;
save('controller','Ak','Bk','Ck','Dk');

Ap = A;
Bp = B;
Cp = C;
Hp = E;
save('plant','Ap','Bp','Cp','Hp');