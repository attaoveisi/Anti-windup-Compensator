function out = opt_par(inp)

load identified_t1
load identified_t1_full

Ap = A;
Bp = B;
Cp = C;
Hp = H;
Dp = zeros(1,2);

np = sqrt(numel(Ap));
nu = numel(Bp)/np;
nw = numel(Hp)/np;
ny = numel(Cp)/np;

load controller
nk = sqrt(numel(Ak));

% inp = [389193.8230150841	516431.97294888785	326177.42254033213	35015.021665047025	248798.90010187373	188268.6168292306	231282.07810577843	67476.80305109014	152835.64225245814	85869.93312043496	326.09913247742935	6912.7639703674295	54930.0094723954	76204.43353884663];
Qp = diag(inp(1,1:np));
Rp = diag(inp(1,np+1:end));

%% LMI
% Define variables
Q = sdpvar(np,np);
U = sdpvar(nu,nu);
X2 = sdpvar(nu,nu,'full');
X1 = sdpvar(nu,np);
gama = sdpvar(1);

% Define constraints and objective
Pi1 = [Q*A'+A*Q,B*U+X1';U*B'+X1,X2'+X2-2*U];
Pi2 = [gama*eye(np),eye(np);eye(np),Q];
Pi3 = [Q*A'+A*Q+B*X1+X1'*B',Q,X1';Q,-Qp,zeros(np,nu);X1,zeros(nu,np),-Rp];

Constraints = [Pi1 <= 0,Pi2 >= 0,Pi3 <= 0,gama >= 0,Q >= 0,U >= 0];
Objective = [gama];

% Set some options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','mosek');

sol = optimize(Constraints,Objective,options);


if sol.problem == 0
    % Extract and display value
    Q = value(Q);
    U = value(U);
    X2 = value(X2);
    X1 = value(X1);

    K = X1*Q^-1;
    L = X2*U^-1;
    options = simset('SrcWorkspace','current');
    sim('AW_test',[],options);
    out = opt_out(end,3)-opt_out(end,2);
else
    
%      display('Hmm, something went wrong!');
     out = 100;
%      yalmiperror(sol.problem)
end
