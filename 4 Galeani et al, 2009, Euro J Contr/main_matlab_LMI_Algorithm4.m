clear 
close all
clc

%% PLANT MODEL
% load identified
% load identified_freq
% global D_hat B_hat P11 Bp C_hat Cy Q11 A_hat Ap
% load identified_t
% B = B/1e4;
load identified_t1
load identified_t1_full

Ap = A;
Bpu = B;
Bpw = H;
Cp = C;
Dpu = zeros(1,2);
Dpw = zeros(1,1);
Cz = C;
Dzu = Dpu;
Dzw = zeros(1,1);

np = sqrt(numel(Ap));
m = numel(Bpu)/np;
q = numel(Bpw)/np;
p = numel(Cp)/np;
l = p;

save plant Ap Bpu Bpw Cp Dpu Dpw Cp Dzw Dzu

%% Controller
%%%%%%% Kalman Filter Design
sysk = ss(Ap,[Bpu Bpw],Cp,0);

[kest,L,P] = kalman(sysk,1e3,1,0);

%%%%%%% Optimal Gain
syst = ss(Ap,Bpu,Cp,0);
[K,S,e] = lqry(syst,1e4,1e0*eye(m),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ac = rlqg.A;
Bc = rlqg.B;
Cc = rlqg.C;
Dc = rlqg.D;
Bcw = zeros(np,q);
Dcw = zeros(m,q);
nc = sqrt(numel(Ac));
save controller.mat Ac Bc Cc Dc Bcw Dcw

%%
Delta = eye(m,m)-Dc*Dpu;
AA = [Ap+Bpu*inv(Delta)*Dc*Cp Bpu*inv(Delta)*Cc;Bc*(eye(p,p)+Dpu*inv(Delta)*Dc)*Cp Ac+Bc*Dpu*inv(Delta)*Cc];
B2 = [Bpu*inv(Delta)*(Dcw+Dc*Dpw)+Bpw;Bc*Dpu*inv(Delta)*(Dcw+Dc*Dpw)+Bcw+Bc*Dpw];
C2 = [Cz+Dzu*inv(Delta)*Dc*Cp Dzu*inv(Delta)*Cc];
D22 = Dzw+Dzu*inv(Delta)*(Dcw+Dc*Dpw);

Bv = [Bpu*inv(Delta)*[zeros(m,nc) eye(m,m)];Bc*Dpu*inv(Delta)*[zeros(m,nc) eye(m,m)]+[eye(nc,nc) zeros(nc,m)]];
Bphi = [Bpu*(eye(m,m)+inv(Delta)*Dc*Dpu);Bc*Dpu*(eye(m,m)+inv(Delta)*Dc*Dpu)];
C1 = [inv(Delta)*Dc*Cp inv(Delta)*Cc];
Cv1 = inv(Delta)*[zeros(m,nc) eye(m,m)];
D1 = inv(Delta)*Dc*Dpu;
Cv2 = Dzu*inv(Delta)*[zeros(m,nc) eye(m,m)];
D2 = Dzu*(eye(m,m)+inv(Delta)*Dc*Dpu);

%% Eq. (24)
% LMI Definition
% Define variables
naw = np+nc;
n = naw+np+nc;

Aaw = zeros(naw,naw);
Baw = zeros(naw,m);
Caw = zeros(nc+m,naw);

Daw = sdpvar(nc+m,m);
gama = sdpvar(1,1);
Q = sdpvar(n,n);

% This should be defined manually
S = eye(m,m);

Ai = [AA Bv*Caw;zeros(naw,np+nc) Aaw];
B1i = [Bphi+Bv*Daw;Baw];
B2i = [B2;zeros(naw,q)];
C1i = [C1 Cv1*Caw];
C2i = [C2 Cv2*Caw];
D21i = D2+Cv2*Daw;
D22i = Dzw+Dzu*inv(Delta)*(Dcw+Dc*Dpw);
D11i = D1+Cv1*Daw;
D12i = inv(Delta)*(Dcw+Dc*Dpw);

% Define constraints and objective
Pi1 = [Ai*Q+Q'*Ai',B1i*S-Q*C1i',B2i,Q*C2i';(B1i*S-Q*C1i')',-S-S'-D11i*S-S'*D11i',-D12i,S*D21i';B2i',-D12i',-eye(q,q),D22i';C2i*Q',D21i*S',D22i,-gama*eye(l,l)];

Constraints2 = [Pi1 <= 0,gama >= 0,Q >= 0];
Objective2 = [];

%% Finding Feasible Solution &&&&& Test
options2 = sdpsettings('verbose',1,'solver','mosek','radius',6e4);
% Solve the problem
sol2 = optimize(Constraints2,Objective2,options2);
% sol2 = optimize(Constraints2);

% Analyze error flags
if sol2.problem == 0
    
    Daw = value(Daw);
    gama = value(gama);

else
    
    display('Something went wrong!');
    sol2.info
    yalmiperror(sol2.problem)
    
end

Khat = ss(Ac,[Bc eye(nc) zeros(nc,m)],Cc,[Dc zeros(m,nc) eye(m)]);