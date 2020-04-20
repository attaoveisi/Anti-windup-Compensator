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
Y = sdpvar(nk,nu);
gama = sdpvar(1);
k = 0;

% Define constraints and objective
Pi1 = [F'*Q1+Q1*F-H'*Y'-Y*H,Q1*G-Y*L,k*Y,H';G'*Q1-L'*Y',-gama*eye(ny,ny),zeros(ny,nu),L';-k*Y',zeros(nu,ny),-gama*eye(nu,nu),k*eye(nu,nu);H,L,k*eye(nu,nu),-gama*eye(nu,nu)];

Constraints = [Pi1 <= 0,gama >= 1000,Q1 >= 0];
Objective = [];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek','radius',1e5);

sol = optimize(Constraints,Objective,options);
if sol.problem == 0
    Q1 = value(Q1);
    Y = value(Y);
    gama = value(gama);
    Pi1 = [F'*Q1+Q1*F-H'*Y'-Y*H,Q1*G-Y*L,k*Y,H';G'*Q1-L'*Y',-gama*eye(ny,ny),zeros(ny,nu),L';-k*Y',zeros(nu,ny),-gama*eye(nu,nu),k*eye(nu,nu);H,L,k*eye(nu,nu),-gama*eye(nu,nu)];

    eig(Q1)
    eig(Pi1)
    gama
else
    sol.info
    error('Sth is wrong')
end

%%
M = inv(Q1)*Y;