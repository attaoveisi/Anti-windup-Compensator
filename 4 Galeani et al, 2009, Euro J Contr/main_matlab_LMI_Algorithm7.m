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
Qp = eye(np,np);
Rp = eye(m,m);

%% LMI Definition
% Define variables
Q = sdpvar(np,np);
X1 = sdpvar(m,np);
X2 = sdpvar(m,m,'full');
gama = sdpvar(1);
U = sdpvar(m,m);

% Define constraints and objective
Pi1 = [Q*Ap'+Ap*Q,Bpu*U+X1';(Bpu*U+X1')',X2'+X2-2*U];
Pi2 = [gama*eye(np,np),eye(np,np);eye(np,np),Q];
Pi3 = [Ap*Q+Q*Ap'+Bpu*X1+X1'*Bpu',Q',X1';Q,-inv(Qp),zeros(np,m);X1,zeros(m,np),-inv(Rp)];

Constraints = [Pi1 <= 0,Pi2 >= 0,Pi3 <= 0,gama >= 0,U >= 0,Q >= 0];

Objective = [gama];

%% Finding Feasible Solution &&&&& Test
options = sdpsettings('verbose',1,'solver','mosek');
% Solve the problem
sol = optimize(Constraints,Objective,options);
% Analyze error flags
if sol.problem == 0
    
    X1 = value(X1);
    X2 = value(X2);
    gama = value(gama);
    U = value(U);
    Q = value(Q);
    
else
    
    display('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
    
end

%%
K = X1*inv(Q);
L = X2*inv(U);

Khat = ss(Ac,[Bc zeros(nc,m)],Cc,[Dc eye(m)]);

ub = 100;