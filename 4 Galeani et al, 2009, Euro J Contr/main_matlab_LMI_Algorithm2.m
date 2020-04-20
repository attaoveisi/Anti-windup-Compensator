clear 
close all
clc

% global Q AA Bv naw np nc Bphi B2 q C1 Cv1 D1 Cv2 D2 C2 Dzw Dzu Delta Dcw Dc Dpw gama Aaw Baw Caw Daw m l
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

%% LMI Definition
% Define variables
X = sdpvar(np+nc,np+nc);
Y = sdpvar(np+nc,np+nc);
gama = sdpvar(1);
Y1 = Y(np,np);

% Define constraints and objective
Pi1 = [AA'*X+X*AA X*B2 C2';B2'*X -eye(q,q) D22';C2 D22 -gama*eye(l,l)];
Pi2 = [Y1*Ap'+Ap*Y1 Bpw Y1*Cz';Bpw' -eye(q,q) Dzw';Cz*Y1 Dzw -gama*eye(l,l)];
Pi3 = [X eye(np+nc,np+nc);eye(np+nc,np+nc) Y];

Constraints = [Pi1 <= 0,Pi2 <= 0,Pi3 >= 0,gama <= 755,X >= 0, Y >= 0];
% Constraints = [Pi1 <= -1e-6,Pi2 <= -1e-6,Pi3 >= 1e-6,gama >= 1e-6,X >= 1e-6, Y >= 1e-6];
Objective = [];

%% Finding Feasible Solution &&&&& Test
options = sdpsettings('verbose',1,'solver','mosek');
% Solve the problem
sol = optimize(Constraints,Objective,options);
% Analyze error flags
if sol.problem == 0
    
    X = value(X);
    Y = value(Y);
    gama = value(gama);
    Y1 = Y(np,np);
    Pi1 = [AA'*X+X*AA X*B2 C2';B2'*X -eye(q,q) D22';C2 D22 -gama*eye(l,l)];
    Pi2 = [Y1*Ap'+Ap*Y1 Bpw Y1*Cz';Bpw' -eye(q,q) Dzw';Cz*Y1 Dzw -gama*eye(l,l)];
    Pi3 = [X eye(np+nc,np+nc);eye(np+nc,np+nc) Y];
    
    real(eig(Pi1))
    real(eig(Pi2))
    real(eig(Pi3))
    real(eig(X))
    real(eig(Y))
    gama

else
    
    display('Hmm, something went wrong!');
    sol.info
    yalmiperror(sol.problem)
    
end

%%
[Mt,N] = lu(eye(np+nc,np+nc)-X*Y);
M = Mt';

test1 = max(max(abs(M'*N-(eye(np+nc)-X*Y))))
Q = [Y,eye(np+nc,np+nc);N,zeros(np+nc,np+nc)]*inv([eye(np+nc,np+nc),X;zeros(np+nc,np+nc),M]);
eig(Q)

%% Eq. (19)
% LMI Definition
% Define variables
naw = np+nc;

Bv = [Bpu*inv(Delta)*[zeros(m,nc) eye(m,m)];Bc*Dpu*inv(Delta)*[zeros(m,nc) eye(m,m)]+[eye(nc,nc) zeros(nc,m)]];
Bphi = [Bpu*(eye(m,m)+inv(Delta)*Dc*Dpu);Bc*Dpu*(eye(m,m)+inv(Delta)*Dc*Dpu)];
C1 = [inv(Delta)*Dc*Cp inv(Delta)*Cc];
Cv1 = inv(Delta)*[zeros(m,nc) eye(m,m)];
D1 = inv(Delta)*Dc*Dpu;
Cv2 = Dzu*inv(Delta)*[zeros(m,nc) eye(m,m)];
D2 = Dzu*(eye(m,m)+inv(Delta)*Dc*Dpu);

S = 2e2*eye(m,m);

Aaw = sdpvar(naw,naw,'full');
Baw = sdpvar(naw,m);
Caw = sdpvar(nc+m,naw);
Daw = sdpvar(nc+m,m);

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
Pi4 = [Ai*Q+Q'*Ai',B1i*S-Q*C1i',B2i,Q*C2i';(B1i*S-Q*C1i')',-S-S'-D11i*S-S'*D11i',-D12i,S*D21i';B2i',-D12i',-eye(q,q),D22i';C2i*Q',D21i*S',D22i,-gama*eye(l,l)];

Constraints2 = [Pi4 <= 0,Aaw <= 0];
Objective2 = [];

% Finding Feasible Solution &&&&& Test
options2 = sdpsettings('verbose',1,'solver','sedumi');
% Solve the problem
sol2 = optimize(Constraints2,Objective2,options2);
% sol2 = optimize(Constraints2);

% Analyze error flags
if sol2.problem == 0
    
    Aaw = value(Aaw);
    Baw = value(Baw);
    Caw = value(Caw);
    Daw = value(Daw);
%     gama = value(gama);

else
    
    display('Something went wrong!');
    sol2.info
    yalmiperror(sol2.problem)
    
end
