clear 
close all
clc

% %% PLANT MODEL
% % load identified
% % load identified_freq
% load identified_t
% B = B/1e4;
% 
% D = zeros(1,2);
% save plant A B H C D
% 
% n = sqrt(numel(A));
% nu = numel(B)/n;
% nw = numel(H)/n;
% ny = numel(C)/n;
% 
% %% Controller
% %%%%%%% Kalman Filter Design
% sysk = ss(A,[B H],C,0);
% 
% [kest,~,~] = kalman(sysk,1e3,1,0);
% 
% %%%%%%% Optimal Gain
% syst = ss(A,B,C,0);
% [K,S,e] = lqry(syst,1e0,1e0*eye(nu),0);
% 
% %%%%%%% Controller
% rlqg = lqgreg(kest,K);
% 
% F = rlqg.A;
% G = rlqg.B;
% H = rlqg.C;
% L = rlqg.D;
% save controller.mat F G H L
% nk = sqrt(numel(F));

%%
G = tf(1,[75 1])*[0.878 -0.864;1.082 -1.096];
Gss = ss(G);
np = sqrt(numel(Gss.A));
A = Gss.A;
B = Gss.B;
C = Gss.C;
D = Gss.D;
nu = numel(B)/np;
ny = numel(C)/np;
K = tf([75 1],[1.43 0])*[45.38 -35.77;44.8 -36.23];

E = 1/75*[0.878 -0.864;1.082 -1.096];
H = E'*E;
T = inv(H);

%% LMI
% Define variables
Q = sdpvar(np,np);
X = sdpvar(nu,np);
gama = sdpvar(1);
u = sdpvar(nu,1);
U = diag(u);  

% Define constraints and objective
Pi1 = [A*Q+Q*A'+B*X+X'*B',B*U*T-X',zeros(np,nu),Q*C'+X'*D';T*U*B'-X,-U*T-T*U,eye(nu,nu),T*U*D';zeros(nu,np),eye(nu,nu),-gama*eye(nu,nu),zeros(nu,ny);C*Q+D*X,D*U*T,zeros(ny,nu),-gama*eye(ny,ny)];

Constraints = [Pi1 <= 0,gama >= 0,Q >= 0,U >= 0];
Objective = [];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

sol = optimize(Constraints,Objective,options);

if sol.problem ~= 0
    sol.info
    yalmiperror(sol.problem)
    error('No answer')
end
    
Q = value(Q);
gama = value(gama);
X = value(X);
u = value(u);
U = diag(u);  

Pi1 = [A*Q+Q*A'+B*X+X'*B',B*U*T-X',zeros(np,nu),Q*C'+X'*D';T*U*B'-X,-U*T-T*U,eye(nu,nu),T*U*D';zeros(nu,np),eye(nu,nu),-gama*eye(nu,nu),zeros(nu,ny);C*Q+D*X,D*U*T,zeros(ny,nu),-gama*eye(ny,ny)];

test = eig(Pi1)
%%
F = X*inv(Q);
