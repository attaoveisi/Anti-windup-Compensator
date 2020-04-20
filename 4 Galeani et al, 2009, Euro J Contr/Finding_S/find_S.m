function output = find_S(input)
global Q AA Bv naw np nc Bphi B2 q C1 Cv1 D1 Cv2 D2 C2 Dzw Dzu Delta Dcw Dc Dpw gama Aaw Baw Caw Daw m l

S = diag(input);

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

Constraints2 = [Pi4 <= 0,gama >= 0];
Objective2 = gama;

% Finding Feasible Solution &&&&& Test
options2 = sdpsettings('verbose',0,'solver','mosek');
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

% else
    
%     display('Something went wrong!');
%     sol2.info
%     yalmiperror(sol2.problem)
    
end
output = sol2.problem;