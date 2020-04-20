clear 
close all
clc

%% PLANT MODEL
load identified4

load SS_model
A = A_ss_retain_modes_model_1;
B = (B_ss_modal_retain_DOFs_retain_modes(:,3:4));
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
lb = -1e0*ones(1,nvars);
ub = 1e0*ones(1,nvars);
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

lmiterm([-3 1 1 Q2],1,1);

lmiterm([-4 1 1 gama_w],1,1);

lmiterm([-5 1 1 gama_d],1,1);

lmiterm([-6 1 1 gama_v],1,1);

LMISYS = getlmis;

% Finding Feasible Solution 
n_dec = decnbr(LMISYS);
c = zeros(1,n_dec);
c(1:end) = 0;
c(1:end-1) = 1;
c(1:end-2) = 0;
options = [1e-15 1000 40 0 0];
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
result_opt = [eig(A) eig(N) eig(Ak) eig(A+B*T)]
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