clear 
close all
clc

%% PLANT MODEL
% load identified
% load identified_freq
global D_hat B_hat P11 Bp C_hat Cy Q11 A_hat Ap
load identified_t1
load identified_t1_full
% load identified_t
% B = B/1e4;

Ap = A;
Bp = B;
Bw = H;
Cy = C;
Dyw = zeros(1,1);
Cz = C;
Dzw = zeros(1,1);
Dz = zeros(1,2);

np = sqrt(numel(Ap));
nu = numel(Bp)/np;
nw = numel(Bw)/np;
ny = numel(Cy)/np;
nz = ny;

save plant Ap Bp Bw Cy Dyw Cz Dzw Dz

%%  MODEL
Li = eye(nz,nz);
Rj = eye(nw,nw);

%% LMI Definition
lambda = 1;

setlmis([]);

Q11 = lmivar(1,[np 1]);
P11 = lmivar(1,[np 1]);
A_hat = lmivar(2,[np np]);
B_hat = lmivar(2,[np ny]);
C_hat = lmivar(2,[nu np]);
D_hat = lmivar(2,[nu ny]);
GAMA1_hat = lmivar(2,[np nu]);
GAMA2_hat = lmivar(2,[nu nu]);
miu = lmivar(1,[1 1]);
ro = lmivar(1,[1 1]);

%LMI terms
lmiterm([1 1 1 Q11],Ap,1,'s');
lmiterm([1 1 1 C_hat],Bp,1,'s');
lmiterm([1 1 1 Q11],lambda,1);
lmiterm([1 1 2 0],Ap+lambda*eye(size(Ap)));
lmiterm([1 1 2 D_hat],Bp,Cy);
lmiterm([1 1 2 -A_hat],1,1);
lmiterm([1 1 3 0],Bw*Rj);
lmiterm([1 1 3 D_hat],Bp,Dyw*Rj);
lmiterm([1 1 4 GAMA2_hat],-Bp,1);
lmiterm([1 1 4 -C_hat],1,1);
lmiterm([1 2 2 P11],1,Ap,'s');
lmiterm([1 2 2 B_hat],1,Cy,'s');
lmiterm([1 2 2 P11],lambda,1);
lmiterm([1 2 3 P11],1,Bw*Rj);
lmiterm([1 2 3 B_hat],1,Dyw*Rj);
lmiterm([1 2 4 -D_hat],Cy',1);
lmiterm([1 2 4 GAMA1_hat],-1,1);
lmiterm([1 3 3 miu],-1,1);
lmiterm([1 3 4 -D_hat],Rj'*Dyw',1);
lmiterm([1 4 4 GAMA2_hat],-1,1,'s');

lmiterm([-2 1 1 Q11],1,1);
lmiterm([-2 1 2 0],1);
lmiterm([-2 2 2 P11],1,1);

lmiterm([-3 1 1 miu],1,1);

lmiterm([-4 1 1 Q11],1,1);

lmiterm([-5 1 1 P11],1,1);

lmiterm([-6 1 1 Q11],lambda,1);
lmiterm([-6 1 2 0],lambda);
lmiterm([-6 1 3 0],0);
lmiterm([-6 1 4 -Q11],1,Cz'*Li');
lmiterm([-6 1 4 -C_hat],1,Dz'*Li');
lmiterm([-6 1 5 GAMA2_hat],Bp,1);
lmiterm([-6 1 5 -C_hat],-1,1);
lmiterm([-6 2 2 P11],lambda,1);
lmiterm([-6 2 3 0],0);
lmiterm([-6 2 4 -Q11],1,Cz'*Li');
lmiterm([-6 2 4 -D_hat],Cy',Dz'*Li');
lmiterm([-6 2 5 GAMA1_hat],1,1);
lmiterm([-6 2 5 -D_hat],Cy',1);
lmiterm([-6 3 3 ro],1,1);
lmiterm([-6 3 3 miu],-1,1);
lmiterm([-6 3 4 0],Rj'*Dzw'*Li');
lmiterm([-6 3 4 -D_hat],Rj'*Dyw',Dzw'*Li');
lmiterm([-6 4 4 ro],1,1);
lmiterm([-6 4 5 GAMA2_hat],-Li*Dz,1);
lmiterm([-6 5 5 GAMA2_hat],1,1,'s');

lmiterm([-7 1 1 ro],1,1);

lmiterm([-8 1 1 ro],1,1);
lmiterm([-8 1 1 miu],-1,1);

lmiterm([9 1 1 ro],1,1);
lmiterm([9 1 1 0],-500);

LMISYS = getlmis;

%% Finding Feasible Solution 
n_dec = decnbr(LMISYS);
c = zeros(1,n_dec);
% c(end-1) = 1;
% c(end) = 1;
% options = [1e-6 1000 0 0 0];
% [~,x] = mincx(LMISYS,c,options,[],[]);
% [~,x] = feasp(LMISYS,[0 1000 -1 99 0]);
[~,x] = feasp(LMISYS);
Q11 = dec2mat(LMISYS,x,1);
P11 = dec2mat(LMISYS,x,2);
A_hat = dec2mat(LMISYS,x,3);
B_hat = dec2mat(LMISYS,x,4);
C_hat = dec2mat(LMISYS,x,5);
D_hat = dec2mat(LMISYS,x,6);
GAMA1_hat = dec2mat(LMISYS,x,7);
GAMA2_hat = dec2mat(LMISYS,x,8);
miu = dec2mat(LMISYS,x,9);
ro = dec2mat(LMISYS,x,10);

%% Method 1
[Q12,P12t] = lu(eye(np,np)-Q11*P11);
P12 = P12t';
test1 = max(max(abs(Q12*P12t-(eye(np,np)-Q11*P11))))
Dk = D_hat;
Bk = inv(P12)*(B_hat-P11*Bp*Dk);
Ck = (C_hat-Dk*Cy*Q11)*inv(Q12');
Ak = inv(P12)*(A_hat-P12*Bk*Cy*Q11-P11*Bp*Ck*Q12'-P11*(Ap+Bp*Dk*Cy)*Q11)*inv(Q12');
nk = sqrt(numel(Ak));

%%
Ai = [Ap+Bp*Dk*Cy Bp*Ck;Bk*Cy Ak];
Bwi = [Bw+Bp*Dk*Dyw;Bk*Dyw];
Bvi = [-Bp;zeros(nk,nu)];
Bksii = [zeros(np,nk) Bp;eye(nk,nk) zeros(nk,nu)];
Cui = [Dk*Cy Ck];
Duwi = Dk*Dyw;
Duksii = [zeros(nu,nk) eye(nu,nu)];
Czi = [Cz+Dz*Dk*Cy Dz*Ck];
Dzwi = Dzw+Dz*Dk*Dyw;
Dzvi = -Dz;
Dzksii = [zeros(nz,nk) Dz];

%%
% Define variables
Q22 = sdpvar(nk,nk);
M = sdpvar(nu,nu,'full');
X = sdpvar(nk+nu,nu);

Q = [Q11,Q12;Q12',Q22];

% Define constraints and objective
Pi1 = [Q*Ai'+Ai*Q,Bvi*M'+Q'*Cui';M*Bvi'+Cui*Q,-2*M]+([-Bksii;-Duksii]*X*[zeros(nu,np+nk) eye(nu,nu)])+([-Bksii;-Duksii]*X*[zeros(nu,np+nk) eye(nu,nu)])';

Constraints = [Pi1 <= 0,M >= 0,Q >= 0];

optimize(Constraints);

%%
Q22 = value(Q22);
M = value(M);
X = value(X);
Q = [Q11,Q12;Q12',Q22];
W = inv(M);

%% Check
Q = [Q11,Q12;Q12',Q22];
Pi1 = [Q*Ai'+Ai*Q,Bvi*M'+Q'*Cui';M*Bvi'+Cui*Q,-2*M]+([-Bksii;-Duksii]*X*[zeros(nu,np+nk) eye(nu,nu)])+([-Bksii;-Duksii]*X*[zeros(nu,np+nk) eye(nu,nu)])';
eig(Q)
eig(Pi1)
eig(M)

%%
GAMA2 = GAMA2_hat*W-eye(size(GAMA2_hat*W));
GAMA1 = inv(P12)*(GAMA1_hat-P11*Bp*(GAMA2+eye(size(GAMA2)))*M)*W;
