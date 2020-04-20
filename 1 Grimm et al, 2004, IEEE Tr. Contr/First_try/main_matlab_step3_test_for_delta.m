clear 
close all
clc

%% Load variables
load variables_for_step3
load plant
load controller

nc = nn;
nu = 2;

for delta = 1e14:1e14:1e17
    % delta=1e16;
    Vs=eye(nu,nu);
    U=delta*(Vs);
    Sai=[Q*Ao'+Ao*Q,Bqo*U+Q*Cyo',Bw,Q*Czo';U*Bqo'+Cyo*Q,Dyqo*U+U*Dyqo'-2*U,Dyw,U*Dzqo';Bw',Dyw',-gama*eye(nw,nw),Dzw';Czo*Q,Dzqo*U,Dzw,-gama*eye(nz,nz)];
    H=[H1,H2,zeros(naw+nv,nw),H3];
    G=[G1*Q,G2*U,zeros(naw+nu,nw+nz)];

    % LMI Definition
    setlmis([]);
    GAMA = lmivar(2,[naw+nv naw+nu]);
    %LMI terms
    lmiterm([1 1 1 GAMA],H',G,'s'); % LMI #1:
    lmiterm([1 1 1 0],Sai); % LMI #1:
    LMISYS = getlmis;
    [~,x] = feasp(LMISYS,[0 0 0 0 1]);
    if x ~= []
        GAMA = dec2mat(LMISYS,x,1);
        GAMA1 = GAMA(1:naw,1:naw);
        GAMA2 = GAMA(1:naw,naw+1:end);
        GAMA3 = GAMA(naw+1:end,1:naw);
        GAMA4 = GAMA(naw+1:end,naw+1:end);

        save GAMA GAMA1 GAMA2 GAMA3 GAMA4
        break
    end
end