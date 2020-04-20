clear 
close all
clc

%% PLANT MODEL
% load identified
% load identified_freq
% load identified_t
% B = B/1e4;
load identified_t1
load identified_t1_full

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = zeros(1,2);
save plant Ap Bp Hp Cp Dp

np = sqrt(numel(Ap));
nu = numel(Bp)/np;
nw = numel(Hp)/np;
ny = numel(Cp)/np;

%% Controller
%%%%%%% Kalman Filter Design
sysk = ss(Ap,[Bp Hp],Cp,0);

[kest,L,P] = kalman(sysk,5e3,1,0);

%%%%%%% Optimal Gain
syst = ss(Ap,Bp,Cp,0);
[K,S,e] = lqry(syst,1e4,1e0*eye(nu),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ak = rlqg.A;
Bk = rlqg.B;
Ck = rlqg.C;
Dk = rlqg.D;
save controller.mat Ak Bk Ck Dk
nk = sqrt(numel(Ak));

%% LMI
% Define variables
Q = sdpvar(np,np);
U = sdpvar(nu,nu);
X2 = sdpvar(nu,nu,'full');
X1 = sdpvar(nu,np);
gama = sdpvar(1);
Qp1 = sdpvar(np,1);
Rp1 = sdpvar(nu,1);
Qp = diag(Qp1);
Rp = diag(Rp1);
% Qp = diag([1,1,1,1,1,1,1,10,10,10,10,1]);
% % QP = C'*C;
% Rp = eye(nu,nu);
% inp = [389193.8230150841 516431.97294888785 326177.42254033213 35015.021665047025 248798.90010187373 188268.6168292306 231282.07810577843 67476.80305109014 152835.64225245814 85869.93312043496 326.09913247742935 6912.7639703674295 54930.0094723954 76204.43353884663];
% Qp = diag(inp(1,1:np));
% Rp = diag(inp(1,np+1:end));
% Qp = 1e2*Cp'*Cp;
% Rp = eye(2);

% Define constraints and objective
Pi1 = [Q*A'+A*Q,B*U+X1';U*B'+X1,X2'+X2-2*U];
Pi2 = [gama*eye(np),eye(np);eye(np),Q];
Pi3 = [Q*A'+A*Q+B*X1+X1'*B',Q,X1';Q,-Qp,zeros(np,nu);X1,zeros(nu,np),-Rp];

Constraints = [Pi1 <= 0,Pi2 >= 0,Pi3 <= 0,gama >= 0,Q >= 0,U >= 0,Q(:) <= 50,U(:) >= 2000];
Objective = [gama];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

sol = optimize(Constraints,Objective,options);


if sol.problem == 0
    % Extract and display value
    Q = value(Q);
    U = value(U);
    X2 = value(X2);
    X1 = value(X1);
    gama = value(gama);
    Qp1 = value(Qp1);
    Rp1 = value(Rp1);
    Qp = diag(Qp1);
    Rp = diag(Rp1);

    Pi1 = [Q*A'+A*Q,B*U+X1';U*B'+X1,X2'+X2-2*U];
    Pi2 = [gama*eye(np),eye(np);eye(np),Q];
    Pi3 = [Q*A'+A*Q+B*X1+X1'*B',Q,X1';Q,-Qp,zeros(np,nu);X1,zeros(nu,np),-Rp];

    eig(Pi1)
    eig(Pi2)
    eig(Pi3)

    K = X1*Q^-1;
    L = X2*U^-1;
    sim('AW_test')
else
    
%      display('Hmm, something went wrong!');
     sol.info
     yalmiperror(sol.problem)
end
