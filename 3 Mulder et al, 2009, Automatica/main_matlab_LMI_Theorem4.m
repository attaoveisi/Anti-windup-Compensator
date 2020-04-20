clear 
close all
clc

%% PLANT MODEL
% load identified
% load identified_freq
global D_hat B_hat P11 Bp C_hat Cy Q11 A_hat Ap
load identified_t1
load identified_t1_full
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
setlmis([]);

Q11 = lmivar(1,[np 1]);
P11 = lmivar(1,[np 1]);
A_hat = lmivar(2,[np np]);
B_hat = lmivar(2,[np ny]);
C_hat = lmivar(2,[nu np]);
D_hat = lmivar(2,[nu ny]);
GAMA1_hat = lmivar(2,[np nu]);
GAMA2_hat = lmivar(2,[nu nu]);
gama = lmivar(1,[1 1]);

%LMI terms
lmiterm([1 1 1 Q11],Ap,1,'s');
lmiterm([1 1 1 C_hat],Bp,1,'s');
lmiterm([1 1 2 0],Ap);
lmiterm([1 1 2 D_hat],Bp,Cy);
lmiterm([1 1 2 -A_hat],1,1);
lmiterm([1 1 3 0],Bw*Rj);
lmiterm([1 1 3 D_hat],Bp,Dyw*Rj);
lmiterm([1 1 4 -Q11],1,Cz'*Li');
lmiterm([1 1 4 -C_hat],1,Dz'*Li');
lmiterm([1 1 5 GAMA2_hat],-Bp,1);
lmiterm([1 1 5 -C_hat],1,1);
lmiterm([1 2 2 P11],1,Ap,'s');
lmiterm([1 2 2 B_hat],1,Cy,'s');
lmiterm([1 2 3 P11],1,Bw*Rj);
lmiterm([1 2 3 B_hat],1,Dyw*Rj);
lmiterm([1 2 4 0],Cz'*Li');
lmiterm([1 2 4 -D_hat],Cy',Dz'*Li');
lmiterm([1 2 5 GAMA1_hat],-1,1);
lmiterm([1 2 5 -D_hat],Cy',1);
lmiterm([1 3 3 gama],-1,1);
lmiterm([1 3 4 0],Rj'*Dzw'*Li');
lmiterm([1 3 4 -D_hat],Rj'*Dyw',Dz'*Li');
lmiterm([1 3 5 -D_hat],Rj'*Dyw',1);
lmiterm([1 4 4 0],-1);
lmiterm([1 4 5 GAMA2_hat],-Li*Dz,1);
lmiterm([1 5 5 GAMA2_hat],-1,1,'s');

lmiterm([-2 1 1 Q11],1,1);
lmiterm([-2 1 2 0],1);
lmiterm([-2 2 2 P11],1,1);

lmiterm([-3 1 1 gama],1,1);

lmiterm([-4 1 1 Q11],1,1);

lmiterm([-5 1 1 P11],1,1);

lmiterm([6 1 1 gama],1,1);
lmiterm([6 1 1 0],-10);

% lmiterm([6 1 1 A_hat],1,1,'s');

LMISYS = getlmis;

%% Finding Feasible Solution 
n_dec = decnbr(LMISYS);
c = zeros(1,n_dec);
c(end) = 1;
options = [1e-6 1000 0 0 0];
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
gama = dec2mat(LMISYS,x,9);

%% Method 1
[Q12,P12t] = lu(eye(np,np)-Q11*P11);
P12 = P12t';
test1 = max(max(abs(Q12*P12t-(eye(np,np)-Q11*P11))))
Dk = D_hat;
Bk = inv(P12)*(B_hat-P11*Bp*Dk);
Ck = (C_hat-Dk*Cy*Q11)*inv(Q12');
Ak = inv(P12)*(A_hat-P12*Bk*Cy*Q11-P11*Bp*Ck*Q12'-P11*(Ap+Bp*Dk*Cy)*Q11)*inv(Q12');
if sum(real(eig(Ak)) > 0) >= 1
%     error('check line 106');
end
nk = sqrt(numel(Ak));

%% Method 2
% nvars = (np^2)*2;
% 
% [Q120,P12t0] = lu(eye(np,np)-Q11*P11);
% P120 = P12t0';
% k = 1;
% Y_0 = zeros(1,nvars);
% for i = 1:np
%     for j = 1:np
%         Y_0(1,k) = Q120(i,j);
%         k = k+1;
%     end
% end
% for i = 1:np
%     for j = 1:np
%         Y_0(1,k) = P120(i,j);
%         k = k+1;
%     end
% end
% InitialPopulation_Data = Y_0;
% lb = -1e0*ones(1,nvars/2);
% ub = 1e0*ones(1,nvars/2);
% lb = [lb -1e6*ones(1,nvars/2)];
% ub = [ub 1e6*ones(1,nvars/2)];
% PopulationSize_Data = 1000;
% EliteCount_Data = round(PopulationSize_Data/100*2);
% if EliteCount_Data == 0 
%     EliteCount_Data = 1;
% end
% CrossoverFraction_Data = 0.95;
% MigrationFraction_Data = 0.9;
% MigrationInterval_Data = 20;
% Generations_Data = 1000;
% StallGenLimit_Data = 100;
% TolFun_Data = 1e-3;
% TolCon_Data = 1e-3;
% [x,fval,exitflag,output,population,score] = optimization_for_Q12P12(nvars,lb,ub,PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,MigrationInterval_Data,MigrationFraction_Data,Generations_Data,StallGenLimit_Data,TolFun_Data,TolCon_Data,InitialPopulation_Data);
% 
% k = 1;
% for i = 1:np
%     for j = 1:np
%         Q12(i,j) = x(1,k);
%         k = k+1;
%     end
% end
% for i = 1:np
%     for j = 1:np
%         P12(i,j) = x(1,k);
%         k = k+1;
%     end
% end
% 
% test1 = max(max(abs(Q12*P12'-(eye(np,np)-Q11*P11))))
% Dk = D_hat;
% Bk = inv(P12)*(B_hat-P11*Bp*Dk);
% Ck = (C_hat-Dk*Cy*Q11)*inv(Q12');
% Ak = inv(P12)*(A_hat-P12*Bk*Cy*Q11-P11*Bp*Ck*Q12'-P11*(Ap+Bp*Dk*Cy)*Q11)*inv(Q12');
% if sum(real(eig(Ak)) > 0) >= 1
%     warning('Decomposition in Eq. (6) did not work correctly');
% end
% nk = sqrt(numel(Ak));

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
M = sdpvar(nu,nu);
X = sdpvar(nk+nu,nu);

Q = [Q11,Q12;Q12',Q22];

% Define constraints and objective
Pi1 = [Q*Ai'+Ai*Q,Bvi*M'+Q'*Cui';M*Bvi'+Cui*Q,-2*M]+([-Bksii;-Duksii]*X*[zeros(nu,np+nk) eye(nu,nu)])+([-Bksii;-Duksii]*X*[zeros(nu,np+nk) eye(nu,nu)])';

Constraints = [Pi1 <= 0,M >= 0,Q >= 0];
Objective = [gama];
options = sdpsettings('verbose',1,'solver','mosek');

optimize(Constraints,Objective,options);

%%
Q22 = value(Q22);
M = value(M);
X = value(X);
Q = [Q11,Q12;Q12',Q22];
W = inv(M);

GAMA2 = GAMA2_hat*W-eye(size(GAMA2_hat*W));
GAMA1 = inv(P12)*(GAMA1_hat-P11*Bp*(GAMA2+eye(size(GAMA2)))*M)*W;

Lambda = [GAMA1;GAMA2];

P = inv(Q);
Ptest = [P11,P12;P12',P(end-np+1:end,end-np+1:end)];

Khat = ss(Ak,[Bk eye(nk) zeros(nk,nu)],Ck,[Dk zeros(nu,nk) eye(nu)]);