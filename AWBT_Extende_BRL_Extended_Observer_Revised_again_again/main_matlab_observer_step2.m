%% PLANT MODEL
clear
clc
load import2

%% LMI Definition
alfa = 1;
beta = 1;

setlmis([]);

X = lmivar(1,[n 1]);
Th = lmivar(2,[nu n]);
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
lmiterm([1 2 2 0],Q2*P*A+(Q2*P*A)'); 
lmiterm([1 2 2 0],-Kh*C-C'*Kh'); 
lmiterm([1 2 2 0],beta); 
lmiterm([1 2 2 0],-Lh*C-C'*Lh'); 
lmiterm([1 2 3 0],Q2*P*G); 
lmiterm([1 2 3 0],-Kh); 
lmiterm([1 2 3 0],-Lh); 
lmiterm([1 2 4 0],0); 
lmiterm([1 2 5 0],Q2*H); 
lmiterm([1 3 3 gama_w],-1,1);
lmiterm([1 3 4 0],0);
lmiterm([1 3 5 0],0);
lmiterm([1 3 6 0],1);
lmiterm([1 4 4 gama_d],-1,1);
lmiterm([1 5 5 gama_v],-1,1);
lmiterm([1 6 6 0],-alfa);

lmiterm([-2 1 1 X],1,1); 
% lmiterm([-2 2 2 X],-1,1); 
% lmiterm([-2 2 2 0],1e1); 

lmiterm([-3 1 1 gama_w],1,1);

lmiterm([-4 1 1 gama_d],1,1);

lmiterm([-5 1 1 gama_v],1,1);

LMISYS = getlmis;

%% Finding Feasible Solution 
% options = [1e-6 1000 1e6 0 0];
% [copt,xopt] = feasp(LMISYS,options);
[copt,xopt] = feasp(LMISYS);

X = dec2mat(LMISYS,xopt,1);
Th = dec2mat(LMISYS,xopt,2);
gama_w = dec2mat(LMISYS,xopt,3);
gama_d = dec2mat(LMISYS,xopt,4);
gama_v = dec2mat(LMISYS,xopt,5);

T = Th*X;

Ak = N+J*T-L*C;
[eig(A+B*T)]
Bk = L+L*C*H-J*T*H;
Ck = T;
Dk = -T*H;
save('controller','Ak','Bk','Ck','Dk');

Ap = A;
Bp = B;
Cp = C;
Hp = E;
save('plant','Ap','Bp','Cp','Hp');