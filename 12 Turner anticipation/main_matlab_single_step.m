clear 
close force
clc

%% PLANT MODEL
load identified_t1
Hp = H;
Ap = A;
Bp = B;
Cp = C;

np = sqrt(numel(A));

m = numel(B)/np;
nw = numel(H)/np;
ny = numel(C)/np;
nz = ny;

Ap = A;
B2 = B;
C1 = C;
B1 = H;
D11 = zeros(nz,nw);
D12 = zeros(nz,m);
C2 = C1;
D21 = D11;
D22 = D12;
save plant Ap B1 B2 C1 D11 D12 C2 D21 D22

%%
sys1 = ss(Ap,B1,C2,0);

%%%%%%% Kalman Filter Design
sysk = ss(Ap,[B2 B1],C2,0);

[kest,L,P] = kalman(sysk,5e3,1,0);

%%%%%%% Optimal Gain
syst = ss(Ap,B2,C2,0);
[K,~,e] = lqry(syst,1e4,1e0*eye(2),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ac = rlqg.A;
nc = sqrt(numel(Ac));
Bcy = rlqg.B;
Cc = rlqg.C;
Dcy = rlqg.D;
Bcw = zeros(nc,nw);
Dcw = zeros(m,nw);

save controller.mat Ac Bcy Cc Dcy Bcw Dcw

%% Full-order anti-windup compensator construction
% choose Ai
Ai = 0.80*eye(m,m); % !!!!!!!!!!!!!!!!!! NOTE1: There is no solution for arbitrary Ai values

% Create the closed-loop system (Appendix A)
Delta = inv(eye(m,m)-Dcy*D22);
Deltat = inv(eye(nz,nz)-D22*Dcy);
Acl = [Ap+B2*Delta*Dcy*C2 B2*Delta*Cc;Bcy*Deltat*C2 Ac+Bcy*Deltat*D22*Cc];
Bbw = [B1+B2*Delta*(Dcw+Dcy*D21);Bcw+Bcy*Deltat*(D21+D22*Dcw)];
Ccl = [Delta*Dcy*C2 Delta*Cc];
Dbw = Delta*(Dcw+Dcy*D21);
Czcl = [C1+D12*Delta*Dcy*C2 D12*Delta*Cc];
Dbzw = D11+D12*Delta*(Dcw+Dcy*D21);

Bb2 = [-B2*Delta;-Bcy*Deltat*D22];
Bg = [zeros(np,nc) B2*Delta;eye(nc,nc) Bcy*Deltat*D22];
Dt2 = -Delta*Dcy*D22;
Dd = [zeros(m,nc) Delta];
Dtz2 = -D12*Delta;
Dg = [zeros(nz,nc) D12*Delta];

% A = [Acl Bg*Gama3;zeros() Gama1];

%%
% Define variables
S11 = sdpvar(np,np);
S12 = sdpvar(np,nc,'full');
S22 = sdpvar(nc,nc);
V11dum = sdpvar(m,1);
V12dum = sdpvar(m,1);
V22dum = sdpvar(m,1);
Ut1dum = sdpvar(m,1);
V11 = diag(V11dum);
V12 = diag(V12dum);
V22 = diag(V22dum);
Ut1 = diag(Ut1dum);
R11 = sdpvar(np,np);
gama = sdpvar(1);
S = [S11,S12;S12',S22];

R22 = sdpvar(nc,nc);
R12 = sdpvar(np,nc,'full');
R = [R11 R12;R12' R22];

% Define constraints and objective
Pi1 = [R11*Ap'+Ap*R11-2*B2*(Ut1+V11)*B2',-B2*(2*V11+V22+Ut1),B1,R11*C1'-2*B2*(Ut1+V11)*D12';(-B2*(2*V11+V22+Ut1))',-2*(V11+V22-Ai*V12),zeros(m,nw),-(2*V11+V22+Ut1)*D12';B1',zeros(nw,m),-gama*eye(nw,nw),D11';(R11*C1'-2*B2*(Ut1+V11)*D12')',(-(2*V11+V22+Ut1)*D12')',D11,-gama*eye(nz,nz)-2*D12*(Ut1+V11)*D12'];
Pi2 = [S*Acl'+Acl*S,Bb2*Ut1+S*Ccl',Bbw,S*Czcl';(Bb2*Ut1+S*Ccl')',-2*V11-Ut1*Dt2'-Dt2*Ut1,Dbw,Ut1*Dtz2';Bbw',Dbw',-gama*eye(nw,nw),Dbzw';(S*Czcl')',(Ut1*Dtz2')',Dbzw,-gama*eye(nz,nz)];
Pi3 = [R11 R12;R12' R22]-S;
Pi4 = Ai*V11-Ut1;
Pi5 = Ai*V12-V22;

Constraints = [Pi1 <= 0,Pi2 <= 0,Pi3 >= 0,Pi4 <= 0,Pi5 <= 0,R11 >= 0,gama >= 0,Ut1 >= 0,V11 >= 0,V12 >= 0,V22 >= 0,S >= eye(size(S)),S11 >= 0, R22 >= 0];
% Objective = [];
Objective = gama;

% Set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');

% Solve the problem
%%%%%%%%%% Method 1 %%%%%%%%%%%
sol = optimize(Constraints,Objective,options);
% Analyze error flags
if sol.problem == 0
    S11 = value(S11);
    S12 = value(S12);
    S22 = value(S22);
    S = [S11,S12;S12',S22];
    R11 = value(R11);
    gama = value(gama);
    V11dum = value(V11dum);
    V12dum = value(V12dum);
    V22dum = value(V22dum);
    Ut1dum = value(Ut1dum);
    V11 = diag(V11dum);
    V12 = diag(V12dum);
    V22 = diag(V22dum);
    Ut1 = diag(Ut1dum);
    R22 = value(R22);
    R12 = value(R12);
else
 sol.info
end

%% Test the results
Pi1 = [R11*Ap'+Ap*R11-2*B2*(Ut1+V11)*B2',-B2*(2*V11+V22+Ut1),B1,R11*C1'-2*B2*(Ut1+V11)*D12';(-B2*(2*V11+V22+Ut1))',-2*(V11+V22-Ai*V12),zeros(m,nw),-(2*V11+V22+Ut1)*D12';B1',zeros(nw,m),-gama*eye(nw,nw),D11';(R11*C1'-2*B2*(Ut1+V11)*D12')',(-(2*V11+V22+Ut1)*D12')',D11,-gama*eye(nz,nz)-2*D12*(Ut1+V11)*D12'];
Pi2 = [S*Acl'+Acl*S,Bb2*Ut1+S*Ccl',Bbw,S*Czcl';(Bb2*Ut1+S*Ccl')',-2*V11+Ut1*Dt2'-Dt2*Ut1,Dbw,Ut1*Dtz2';Bbw',Dbw',-gama*eye(nw,nw),Dbzw';(S*Czcl')',(Ut1*Dtz2')',Dbzw,-gama*eye(nz,nz)];
Pi3 = R11-S11;
Pi4 = Ai*V11-Ut1;
Pi5 = Ai*V12-V22;

disp('************************************')
disp('Pi1(-)=')
disp(eig(Pi1))
disp('************************************')
disp('Pi2(-)=')
disp(eig(Pi2))
disp('************************************')
disp('Pi3(+)=')
disp(eig(Pi3))
disp('************************************')
disp('Pi4(-)=')
disp(eig(Pi4))
disp('************************************')
disp('Pi5(-)=')
disp(eig(Pi5))
disp('************************************')
disp('R11(+)=')
disp(eig(R11))
disp('************************************')
disp('S(+)=')
disp(eig(S))
disp('************************************')
disp('V11(D+)=')
disp(V11)
disp('************************************')
disp('V12(D+)=')
disp(V12)
disp('************************************')
disp('V22(D+)=')
disp(V22)
disp('************************************')
disp('Ut1(D+)=')
disp(Ut1)
disp('************************************')
disp('gama(+)=')
disp(gama)
disp('************************************')
disp('R-S(+)=')
disp(eig([R11 R12;R12' R22]-S))
disp('************************************')
disp('s^-1PS^-1-S^-1(+)=')
% R12 = S12;
% R22 = S22;
R = [R11,R12;R12',R22];
dumInvS = inv(S);
eig(dumInvS*R*dumInvS-dumInvS)
disp('************************************')

%% Construct R, above Eq. (59)
% Construct Eq. (42)
dumInvS = inv(S);
expression1 = dumInvS*R*dumInvS-dumInvS; % Eq. (59) % !!!!!!!!!!!!!! NOTE2: This expression isn't PD and needs some additional considerations as below
Ps = chol(expression1)';

% Eq. (60)
naw = np+nc;
P33 = eye(naw,naw)+Ps'*S*Ps;

% Eq. (61) - (65)
W11 = V11*((Ut1)^(-2));
W21 = Ut1^(-1)-Ai*W11;
Ut2 = zeros(m,m);
for i = 1:m
    Ut2(i,i) = (V22(i,i)+sqrt((V22(i,i))^2+4*Ai(i,i)*V11(i,i)*V12(i,i)))/(2*Ai(i,i)*Ut1(i,i)*W11(i,i));
end
W12 = V12*Ut2^(-2);
W22 = Ut2^(-1)-Ai*W12;

%% Second LMI in Eq. (35)
% Eq. (A3)
A0 = blkdiag(Acl,zeros(naw,naw));
C0 = [Ccl,zeros(m,naw)];
Cz0 = [Czcl,zeros(nz,naw)];
G2 = [zeros(naw,m);eye(m,m)];
G1 = zeros(naw+m,np+nc+naw);
G1(1:naw,end-naw+1:end) = eye(naw,naw);
H1 = [zeros(naw,np+nc),eye(naw,naw);Bg',zeros(nc+m,naw)];
H2 = [zeros(naw,m);Dd'];
H3 = [zeros(naw,nz);Dg'];

% Define variables
GAMA1 = sdpvar(naw,naw);
GAMA2 = sdpvar(naw,m);
GAMA3 = sdpvar(nc+m,naw);
GAMA4 = sdpvar(nc+m,m);
% GAMA4 = zeros(nc+m,m);
GAMA = [GAMA1,GAMA2;GAMA3,GAMA4];
% gama1 = sdpvar(1);
gama1 = gama;

% Eq. (42)
P = [inv(S) Ps;Ps' P33];

% Eq. (A1): in appendix
A = [Acl Bg*GAMA3;zeros(naw,np+nc) GAMA1];
Bw = [Bbw;zeros(naw,nw)];
Bt1 = [Bg*GAMA4;GAMA2];
Bt2 = [Bb2;zeros(naw,m)];
C = [Ccl Dd*GAMA3];
Dw = Dbw;
Dt1 = Dd*GAMA4;
Cz = [Czcl Dg*GAMA3];
Dzw = Dbzw;
Dtz1 = Dg*GAMA4;

% Eq. (39) & (40)
Wt1 = Ai*W11+W21;
Wt2 = Ai*W12+W22;

H = [H1*P H2*Wt2 H2*Wt1 zeros(nc+m+naw,nw) H3];
G = [G1 G2 zeros(naw+m,m+nw+nz)];
Sai0 = [P*A0 P*Bt2 -P*Bt2 P*Bw zeros(naw+np+nc,nz);
        Wt2*C0 -W22+Wt2*Dt2 -W12-Wt2*Dt2 Wt2*Dw zeros(m,nz);
        Wt1*C0 -W21+Wt1*Dt2 -2*W11-Wt1*Dt2 Wt1*Dw zeros(m,nz);
        zeros(nw,(np+nc+naw)+(m)+(m)) -gama1*eye(nw,nw) zeros(nw,nz);
        Cz0 Dtz2 -Dtz2 Dzw -gama1*eye(nz,nz)];

% Define constraints and objective
% Constraints = [Sai0+Sai0'+H'*GAMA*G+G'*GAMA'*H <= 0,gama1 >= 0];
Constraints = [Sai0+Sai0'+H'*GAMA*G+G'*GAMA'*H <= 0,GAMA1 <= 0];

% Objective = gama1;
Objective = [];

% Solve the problem
%%%%%%%%%% Method 1 %%%%%%%%%%%
sol1 = optimize(Constraints,Objective,options);
% Analyze error flags
if sol1.problem == 0
    GAMA1 = value(GAMA1);
    GAMA2 = value(GAMA2);
    GAMA3 = value(GAMA3);
    GAMA4 = value(GAMA4);
    GAMA = [GAMA1,GAMA2;GAMA3,GAMA4];
%     gama1 = value(gama1);
    
    disp('************************************')
    disp('Constraint(-)=')
    disp(eig(Sai0+Sai0'+H'*GAMA*G+G'*GAMA'*H))
    disp('************************************')
else
 sol1.info
end

%%%%%%%%%% Method 1 %%%%%%%%%%%
% LMI Definition
setlmis([]);

GAMA = lmivar(2,[naw+nc+m naw+m]);

%LMI terms
lmiterm([1 1 1 GAMA],H',G,'s'); % LMI #1:
lmiterm([1 1 1 0],Sai0+Sai0'); % LMI #1:

LMISYS = getlmis;
options = [0 120 1e4 90 0];
[~,x] = feasp(LMISYS,options);
GAMA = dec2mat(LMISYS,x,1);
GAMA1 = GAMA(1:naw,1:naw);
GAMA2 = GAMA(1:naw,naw+1:end);
GAMA3 = GAMA(naw+1:end,1:naw);
GAMA4 = GAMA(naw+1:end,naw+1:end);
GAMA = [GAMA1,GAMA2;GAMA3,GAMA4];
disp(eig(Sai0+Sai0'+H'*GAMA*G+G'*GAMA'*H))
disp('************************************')
disp(eig(GAMA1))