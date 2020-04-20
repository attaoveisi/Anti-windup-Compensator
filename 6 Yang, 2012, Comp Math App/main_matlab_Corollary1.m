clear 
close all
clc

%% PLANT MODEL
% load identified
% load identified_freq
load identified_t1
load identified_t1_full
% load identified_t
% B = B/1e4;

D = zeros(1,2);
save plant A B H C D

n = sqrt(numel(A));
nu = numel(B)/n;
nw = numel(H)/n;
ny = numel(C)/n;

%% Controller
%%%%%%% Kalman Filter Design
sysk = ss(A,[B H],C,0);

[kest,~,~] = kalman(sysk,1e3,1,0);

%%%%%%% Optimal Gain
syst = ss(A,B,C,0);
[K,S,e] = lqry(syst,1e4,1e0*eye(nu),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

F = rlqg.A;
G = rlqg.B;
H = rlqg.C;
L = rlqg.D;
save controller.mat F G H L
nk = sqrt(numel(F));

%% LMI
% Define variables
Q1 = sdpvar(nk,nk);
Q2 = sdpvar(n,n);
Y = sdpvar(nk,nu);
delta1 = sdpvar(1);
m = nu;
x = sdpvar(m,1);
W = diag(x);    % Diagonal matrix

gama = sdpvar(1);
k = 0;

% Define constraints and objective
Z11 = [F'*Q1+Q1*F-H'*Y'-Y*H,Q1*G*C-Y*L*C,zeros(nk,nk);C'*G'*Q1-C'*L'*Y',A'*Q2+Q2*A,zeros(n,nk);zeros(nk,nk),zeros(nk,n),F'*Q1+Q1*F-H'*Y'-Y*H];
Z12 = [Q1*G*D-Y*L*D-H'*W;Q2*B-C'*L'*W;-Y-H'*W];
Z22 = delta1*eye(nu,nu)-2*W-W*L*D-D'*L'*W;

Z = [Z11,Z12;Z12',Z22];

Pi1 = [F'*Q1+Q1*F-H'*Y'-Y*H,Q1*G-Y*L,k*Y,H';G'*Q1-L'*Y',-gama*eye(ny,ny),zeros(ny,nu),L';-k*Y',zeros(nu,ny),-gama*eye(nu,nu),k*eye(nu,nu);H,L,k*eye(nu,nu),-gama*eye(nu,nu)];

Constraints = [Pi1 <= 0,gama >= 0,Z <= 0,delta1 >= 0,W >= 0,Q1 >= 0,Q2 >= 0];
Objective = [];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','sedumi');

sol = optimize(Constraints,Objective,options);
if sol.problem ~= 0
    sol.info
    yalmiperror(sol.problem)
    error('No answer')
end

% optimize(Constraints,Objective);

Q1 = value(Q1);
Y = value(Y);
Q2 = value(Q2);
x = value(x);
W = diag(x); 
delta1 = value(delta1);
gama = value(gama);

Z11 = [F'*Q1+Q1*F-H'*Y'-Y*H,Q1*G*C-Y*L*C,zeros(nk,nk);C'*G'*Q1-C'*L'*Y',A'*Q2+Q2*A,zeros(n,nk);zeros(nk,nk),zeros(nk,n),F'*Q1+Q1*F-H'*Y'-Y*H];
Z12 = [Q1*G*D-Y*L*D-H'*W;Q2*B-C'*L'*W;-Y-H'*W];
Z22 = delta1*eye(nu,nu)-2*W-W*L*D-D'*L'*W;
Z = [Z11,Z12;Z12',Z22];

Pi1 = [F'*Q1+Q1*F-H'*Y'-Y*H,Q1*G-Y*L,k*Y,H';G'*Q1-L'*Y',-gama*eye(ny,ny),zeros(ny,nu),L';-k*Y',zeros(nu,ny),-gama*eye(nu,nu),k*eye(nu,nu);H,L,k*eye(nu,nu),-gama*eye(nu,nu)];

test1 = eig(Z)
test2 = eig(Pi1)
%%
M = inv(Q1)*Y;
