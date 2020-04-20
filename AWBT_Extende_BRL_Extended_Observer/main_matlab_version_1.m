clear 
close all
clc

%% PLANT MODEL
load identified

dum1=size(A);
n=dum1(1,1);
dum2=size(B);
nu=dum2(1,2);
dum3=size(C);
ny=dum3(1,1);
dum4=size(E);
nd=dum4(1,2);

Y=zeros(n,ny);

CEp=((C*E)'*(C*E))^-1;
H=-E*CEp+Y*(eye(ny,ny)-(C*E)*CEp);
P=eye(n,n)+H*C;
J=(eye(n,n)+(-E*CEp+Y*(eye(ny,ny)-(C*E)*CEp))*C)*B;

%% LMI Definition
delta = 0.0001;

setlmis([]);

X = lmivar(1,[n 1]);
Q2 = lmivar(1,[n 1]);
Th = lmivar(2,[nu n]);
Kh = lmivar(2,[n ny]);
gama_w = lmivar(1,[1 1]);
gama_d = lmivar(1,[1 1]);
gama_v = lmivar(1,[1 1]);

%LMI terms
lmiterm([1 1 1 Th],B,2,'s'); 
lmiterm([1 1 1 X],A,2,'s'); 
lmiterm([1 1 2 0],-1);
lmiterm([1 1 2 X],-delta,A');
lmiterm([1 1 2 -Th],-delta,B');
lmiterm([1 1 3 0],0);
lmiterm([1 1 4 0],2*G); 
lmiterm([1 1 5 0],2*E);
lmiterm([1 1 6 0],0);
lmiterm([1 1 7 X],1,C');
lmiterm([1 2 2 0],2*delta);
lmiterm([1 2 3 0],0); % LMI #1: Check -delta*B*T ************ !!!!!!!!! **************
lmiterm([1 2 4 0],-delta*G); 
lmiterm([1 2 5 0],-delta*E); 
lmiterm([1 3 3 Q2],1,P*A,'s'); 
lmiterm([1 3 3 Kh],-1,C,'s'); 
lmiterm([1 3 4 Q2],1,P*G); 
lmiterm([1 3 4 Kh],-1,1); 
lmiterm([1 3 5 0],0); 
lmiterm([1 3 6 Q2],1,H); 
lmiterm([1 4 4 gama_w],-1,1);
lmiterm([1 4 5 0],0);
lmiterm([1 4 6 0],0);
lmiterm([1 4 7 0],1);
lmiterm([1 5 5 gama_d],-1,1);
lmiterm([1 6 6 gama_v],-1,1);
lmiterm([1 7 7 0],-1);

lmiterm([-2 1 1 X],1,1); 

lmiterm([-3 1 1 Q2],1,1);

lmiterm([-4 1 1 gama_w],1,1);

lmiterm([-5 1 1 gama_d],1,1);

lmiterm([-6 1 1 gama_v],1,1);

LMISYS = getlmis;

%% Finding Feasible Solution 
% [~,x] = feasp(LMISYS,[0 1000 -1 99 0]);
[~,x] = feasp(LMISYS);

X = dec2mat(LMISYS,x,1);
Q2 = dec2mat(LMISYS,x,2);
Th = dec2mat(LMISYS,x,3);
Kh = dec2mat(LMISYS,x,4);
gama_w = dec2mat(LMISYS,x,5);
gama_d = dec2mat(LMISYS,x,6);
gama_v = dec2mat(LMISYS,x,7);

K = Q2\Kh;
T = Th/X;
N = P*A-K*C;

[eig(A) eig(A+B*T) eig(N)]

%%
% n_dec = decnbr(LMISYS);
% c = zeros(1,n_dec);
% c(1:end) = 1e0;
% c(1:end-1) = 0;
% c(1:end-2) = 0;
% options = [1e-6 1000 0 0 0];
% [copt,xopt] = mincx(LMISYS,c,options,[],[]);


