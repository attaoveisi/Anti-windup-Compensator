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

savematfile("variables_for_step2_prt2.mat","R","S","naw","np","nc","-v7.3");
//[N,Ndum]=fullrf(R*inv(S)*R-R);
