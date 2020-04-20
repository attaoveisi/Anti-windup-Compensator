clc
clear

loadmatfile('plant.mat');
loadmatfile('controller.mat');
loadmatfile('closed_loop.mat');

Ap = Ap(1:nn,1:nn);
Bpu = Bp(1:nn,:);
Cpy = Cp(:,1:nn);
Cpz = [0,0,1,0];
Bpw = Hp(1:nn,:);

dum1=size(Ap);
np=dum1(1,1);
dum2=size(Bpu);
nu=dum2(1,2);
dum3=size(Cpy);
ny=dum3(1,1);
nz=1;
dum4=size(Bpw);
nw=dum4(1,2);

dum5=size(Ak);
nc=dum5(1,2);

naw=np;
//naw=np+nc;
ncl=np+nc;

Dpyu = zeros(ny,nu);
Dpyw = zeros(ny,nw);
Dpzu = zeros(nz,nu);
Dpzw = zeros(nz,nw);

function [LME, LMI, OBJ]=AW_feas(XLIST)
[R11,R12,R22,S,gama]= XLIST(:)
LME=list(R11-R11',R22-R22',S-S')
LMI=list(-([R11*Ap'+Ap*R11,Bpw,R11*Cpz';Bpw',-gama*eye(nw,nw),Dpzw';Cpz*R11,Dpzw,-gama*eye(nz,nz)]),-([S*Acl'+Acl*S,Bclw,S*Cclz';Bclw',-gama*eye(nw,nw),Dclzw';Cclz*S,Dclzw,-gama*eye(nz,nz)]),([R11,R12;R12',R22]),([R11,R12;R12',R22])-S,S)
OBJ=[gama]
endfunction

R110=eye(np,np);
R120=zeros(np,nc);
R220=zeros(nc,nc);
S0=eye(np+nc,np+nc);
gama0=1e1;

Init_guess=list(R110,R120,R220,S0,gama0);

Mbound=1e5;
abstol=1e-20;
nuu=2;
maxiters=5000;
reltol=1e-20;
options=[Mbound,abstol,nuu,maxiters,reltol];

//Ans_LMI1=lmisolver(Init_guess,AW_feas,options);
Ans_LMI1=lmisolver(Init_guess,AW_feas);

R11=Ans_LMI1(1);
R12=Ans_LMI1(2);
R22=Ans_LMI1(3);
S=Ans_LMI1(4);
gama=Ans_LMI1(5);

R=[R11,R12;R12',R22];

[N,Ndum]=fullrf(R*inv(S)*R-R);

M=eye(naw,naw)+N'*R^-1*N;
Q=[R,N;N',M];

n=naw+nc+np;
nv=nc+nu;

Ao=sysdiag(Acl,zeros(naw,naw));
Bqo=[Bclq;zeros(naw,nu)];
Cyo=[Ccly,zeros(nu,naw)];
Dyqo=Dclyq;
Czo=[Cclz,zeros(nz,naw)];
Dzqo=Dclzq;
H1=[zeros(np+nc,naw),Bclv;eye(naw,naw),zeros(naw,nc+nu)]';
G1=zeros(naw+nu,n);
G1(1:naw,$-naw+1:$)=eye(naw,naw);
G2=[zeros(naw,nu);eye(nu,nu)];
H2=[zeros(nu,naw),Dclyv]';
H3=[zeros(nz,naw),Dclzv]';
Bw=[Bclw;zeros(naw,nw)];
Dzw=Dclzw;
Dyw=Dclyw;

delta=1;
Vs=eye(nu,nu);
U=delta*inv(Vs);
Sai=[Q*Ao'+Ao*Q,Bqo*U+Q*Cyo',Bw,Q*Czo';U*Bqo'+Cyo*Q,Dyqo*U+U*Dyqo'-2*U,Dyw,U*Dzqo';Bw',Dyw',-gama*eye(nw,nw),Dzw';Czo*Q,Dzqo*U,Dzw,-gama*eye(nz,nz)];
H=[H1,H2,zeros(naw+nv,nw),H3];
G=[G1*Q,G2*U,zeros(naw+nu,nw+nz)];

savematfile("variables_for_step3.mat", "H", "G", "Sai", "naw", "nv", "nu", "-v7.3");

function [LME, LMI, OBJ]=AW_feas1(XLIST)
[GAMA1,GAMA2,GAMA3,GAMA4]= XLIST(:)
GAMA=[GAMA1,GAMA2;GAMA3,GAMA4];
LME=[]
LMI=list(-(Sai+G'*GAMA'*H+H'*GAMA*G),-GAMA1)
OBJ=[]
endfunction

GAMA10=[-65.02,198.43,98.11,-66.75;223.94,-697.09,-347.39,247.24;41.17,-98.1,-47.56,55.25;-121.39,309.97,138.31,-131.52];
GAMA20=[0.0688,-0.262,-0.0637,0.1559]';
GAMA30=[41.22,-160.42,-106.41,82.03;-3469.09,8318.57,3423.87,-2388.49;-162.51,386.26,-35.56,71.07;-4584.37,9490.06,-11350.16,11407.08;587.11,-1687.16,-821.25,632.86];
GAMA40=[-0.0622,2.907,0.2338,5.5623,0]';
Init_guess=list(GAMA10,GAMA20,GAMA30,GAMA40);

Ans_LMI2=lmisolver(Init_guess,AW_feas1);

GAMA1=Ans_LMI2(1);
GAMA2=Ans_LMI2(2);
GAMA3=Ans_LMI2(3);
GAMA4=Ans_LMI2(4);

GAMA=[GAMA1,GAMA2;GAMA3,GAMA4];
