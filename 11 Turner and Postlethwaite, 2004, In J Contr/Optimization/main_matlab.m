clear 
close all
clc

%% PLANT MODEL
load ID_results
load identified_250

%%
global Ap Bp Hp Cp Dp Ak Bk Ck Dk Af Bf Hf Cf Df theta
A = sysc_250.A;
H = sysc_250.B(:,4);
B = sysc_250.B(:,1:2);
C = sysc_250.C;
D = sysc_250.D(:,1:2);

Af = sysc_full.A;
Hf = sysc_full.B(:,4);
Bf = sysc_full.B(:,1:2);
Cf = sysc_full.C;
Df = sysc_full.D(:,1:2);

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = D;
% Dp = zeros(1,2);
save plant Ap Bp Hp Cp Dp

np = sqrt(numel(Ap));
nu = numel(Bp)/np;
nw = numel(Hp)/np;
ny = numel(Cp)/np;

%% Controller
%%%%%%% Kalman Filter Design
sysk = ss(Ap,[Bp Hp],Cp,0);

[kest,L,P] = kalman(sysk,5e3,eye(ny),0);

%%%%%%% Optimal Gain
syst = ss(Ap,Bp,Cp,0);
[K,S,e] = lqry(syst,3e2*eye(ny),1e0*eye(nu),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ak = rlqg.A;
Bk = rlqg.B;
Ck = rlqg.C;
Dk = rlqg.D;
save controller.mat Ak Bk Ck Dk
nk = sqrt(numel(Ak));

%%
Delta = (eye(ny)-Dp*Dk)^-1;
Deltat = (eye(nu)-Dk*Dp)^-1;
Ab = [Ap+Bp*Deltat*Dk*Cp Bp*Deltat*Ck;Bk*Delta*Cp Ak+Bk*Delta*Dp*Ck];
B0 = [Bp*Deltat;Bk*Delta*Dp];
Bb = [Bp*Deltat -Bp*Deltat*Dk;Bk*Delta*Dp -Bk*Delta];
Cb1 = [Deltat*Dk*Cp Deltat*Ck];
D01 = Deltat*Dk*Dp;
Db1 = [Deltat -Deltat*Dk];
Cb2 = [Delta*Cp Delta*Dp*Ck];
D02 = Delta*Dp;
Db2 = [Delta*Dp -Delta*Dp*Dk];

%% LMI
% Define variables
Q = sdpvar(np+nk,np+nk);
L = sdpvar(nu+ny,nu);
U_dum = sdpvar(nu,1);
U = diag(U_dum);
gama = sdpvar(1);

% Define constraints and objective
Pi1 = [Q*Ab'+Ab*Q B0*U+Bb*L-Q*Cb1' zeros(np+nk,nu) Q*Cb2';(B0*U+Bb*L-Q*Cb1')' -2*U-U*D01'-D01*U-L'*Db1'-Db1*L,eye(nu,nu) U*D02'+L'*Db2';zeros(nu,np+nk) eye(nu,nu) -gama*eye(nu,nu) zeros(nu,ny);(Q*Cb2')' (U*D02'+L'*Db2')' zeros(ny,nu) -gama*eye(ny,ny)];

% Constraints = [Pi1 <= -10*eye(size(Pi1)),gama >= 0,Q >= 0,U_dum(:) >= 1e4];
Constraints = [Pi1 <= -1e-3*eye(size(Pi1)),gama >= 0,Q >= 0,U_dum(:) >= 0];
Objective = [gama];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

sol = optimize(Constraints,Objective,options);


if sol.problem == 0
    % Extract and display value
    Q = value(Q);
    U_dum = value(U_dum);
    U = diag(U_dum);
    gama = value(gama);
    L = value(L);

    Pi1 = [Q*Ab'+Ab*Q B0*U+Bb*L-Q*Cb1' zeros(np+nk,nu) Q*Cb2';(B0*U+Bb*L-Q*Cb1')' -2*U-U*D01'-D01*U-L'*Db1'-Db1*L,eye(nu,nu) U*D02'+L'*Db2';zeros(nu,np+nk) eye(nu,nu) -gama*eye(nu,nu) zeros(nu,ny);(Q*Cb2')' (U*D02'+L'*Db2')' zeros(ny,nu) -gama*eye(ny,ny)];
    eig(Pi1)
    
    theta = L*inv(U);
%     F1 = tf(400*2*pi,[1 800*2*pi])*eye(nu);
%     F2 = tf(800*2*pi,[1 800*2*pi])*eye(ny);
%     F2 = 1;
%     F1 = 1;

%     sim('AW_test')
else
    
%      display('Hmm, something went wrong!');
     sol.info
%      yalmiperror(sol.problem)
end

