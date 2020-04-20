clear 
close all
clc

%% PLANT MODEL
load identified_t1
load identified_t1_full

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = zeros(1,2);
np = sqrt(numel(Ap));
nu = numel(Bp)/np;
ny = numel(Cp)/np;
save plant Ap Bp Hp Cp Dp

%%
%%%%%%% Kalman Filter Design
sysk=ss(A,[Bp Hp],C,0);

[kest,L,P] = kalman(sysk,5e3,1,0);

%%%%%%% Optimal Gain
syst=ss(A,B,C,0);
[K,S,e] = lqry(syst,1e4,1e0*eye(nu),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ak = rlqg.A;
Bk = rlqg.B;
Ck = rlqg.C;
Dk = rlqg.D;
save controller_reduced.mat Ak Bk Ck Dk

%% PLANT MODEL
load plant
load controller_reduced

nk = sqrt(numel(Ak));
D = Dp;

%%
% SSS = stepinfo(rlqg);
E = dcgain(ss(A,B,C,Dp))';
% E = evalfr(ss(Ak,Bk,Ck,Dk),800).^-1;
H = E'*E;
T = H^-1;

%% LMI
% Define variables
Q = sdpvar(np,np);
X = sdpvar(nu,np);
gama = sdpvar(1);
u = sdpvar(nu,1);
U = diag(u);  

% Define constraints and objective
Pi1 = [A*Q+Q*A'+B*X+X'*B',B*U*T-X',zeros(np,nu),Q*C'+X'*D';T*U*B'-X,-U*T-T*U,eye(nu,nu),T*U*D';zeros(nu,np),eye(nu,nu),-gama*eye(nu,nu),zeros(nu,ny);C*Q+D*X,D*U*T,zeros(ny,nu),-gama*eye(ny,ny)];

Constraints = [Pi1 <= -0.1*eye(size(Pi1)),gama >= 0,Q >= 0,u(:) >= 0.0002];
Objective = [gama];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

sol = optimize(Constraints,Objective,options);

if sol.problem ~= 0
    sol.info
    yalmiperror(sol.problem)
    error('No answer')
end
    
Q = value(Q);
gama = value(gama)
X = value(X);
u = value(u);
U = diag(u);  

Pi1 = [A*Q+Q*A'+B*X+X'*B',B*U*T-X',zeros(np,nu),Q*C'+X'*D';T*U*B'-X,-U*T-T*U,eye(nu,nu),T*U*D';zeros(nu,np),eye(nu,nu),-gama*eye(nu,nu),zeros(nu,ny);C*Q+D*X,D*U*T,zeros(ny,nu),-gama*eye(ny,ny)];

test = eig(Pi1)
test = eig(Q)
%%
F = X*inv(Q);
save matlab1 F E
[eig(A) eig(A+B*F)]

% Mt = ss(A+B*F,B*E.^-1,F,E.^-1)-E.^-1;
% Nt = ss(A+B*F,B*E.^-1,C+Dp*F,Dp*E.^-1);

Mt = ss(A+B*F,B,F,eye(nu))-eye(nu);
Nt = ss(A+B*F,B,C+Dp*F,Dp);

sim('AW_test')