clc
clear

loadmatfile('identified.mat');

dum1=size(A);
n=dum1(1,1);
dum2=size(B);
nu=dum2(1,2);
dum3=size(C);
ny=dum3(1,1);
dum4=size(E);
nd=dum4(1,2);

Y=zeros(n,ny);

CEp=((C*E)'*(C*E))^-1;
H=-E*CEp+Y*(eye(ny,ny)-(C*E)*CEp);
P=eye(n,n)+H*C;
J=(eye(n,n)+(-E*CEp+Y*(eye(ny,ny)-(C*E)*CEp))*C)*B;

function [LME, LMI, OBJ]=AWBT_EBRL(XLIST)
[Q1,Q2,Th,Kh,Qh1,gama_w,gama_d,gama_v]= XLIST(:)
LME=list(Q1-Q1',Q2-Q2',Q1*B-B*Qh1)
LMI=list(-([2*(B*Th+Th'*B'+A'*Q1+Q1*A),zeros(n,n),2*Q1*G,2*Q1*E,zeros(n,ny),C';zeros(n,n),A'*P'*Q2+Q2*P*A-Kh*C-C'*Kh',Q2*P*G-Kh,zeros(n,nd),Q2*H,zeros(n,ny);2*G'*Q1,G'*P'*Q2-Kh',gama_w*eye(ny,ny),zeros(ny,nd+ny),eye(ny,ny);2*E'*Q1,zeros(nd,n+ny),-gama_d*eye(nd,nd),zeros(nd,ny+ny);zeros(ny,n),H'*Q2,zeros(ny,ny+nd),gama_v*eye(ny,ny),zeros(ny,ny);C,zeros(ny,n),eye(ny,ny),zeros(ny,nd+ny),-eye(ny,ny)]),gama_w,gama_d,gama_v,Q1,Q2)
OBJ=[]
endfunction

Q10=1e3*eye(n,n);
Q20=1e0*eye(n,n);
Th0=1e3*rand(nu,n);
Kh0=1e3*rand(n,ny);
Qh10=eye(nu,nu);
gama_w0=1e0;
gama_d0=1e0;
gama_v0=1e0;

Init_guess=list(Q10,Q20,Th0,Kh0,Qh10,gama_w0,gama_d0,gama_v0);

Mbound=1e10;
abstol=5e-6;
nu=2;
maxiters=5000;
reltol=1e-10;
options=[Mbound,abstol,nu,maxiters,reltol];

//Ans_LMI=lmisolver(Init_guess,AWBT_EBRL,options);
Ans_LMI1=lmisolver(Init_guess,AWBT_EBRL);

Q1=Ans_LMI1(1);
Q2=Ans_LMI1(2);
Th=Ans_LMI1(3);
Kh=Ans_LMI1(4);
Qh1=Ans_LMI1(5);
gama_w=Ans_LMI1(6);
gama_d=Ans_LMI1(7);
gama_v=Ans_LMI1(8);

