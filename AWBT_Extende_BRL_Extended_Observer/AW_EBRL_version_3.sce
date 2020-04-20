clc
clear

loadmatfile('identified2.mat');

A=A(1:4,1:4);
B=B(1:4,:);
E=E(1:4,:);
G=G(1:4,:);
C=C(:,1:4);

dum1=size(A);
n=dum1(1,1);
dum2=size(B);
nu=dum2(1,2);
dum3=size(C);
ny=dum3(1,1);
dum4=size(E);
nd=dum4(1,2);

CEp=((C*E)'*(C*E))^-1;
H1=-E*CEp;
H2=eye(ny,ny)-(C*E)*CEp;
P1=eye(n,n)-E*CEp*C;
P2=(eye(ny,ny)-(C*E)*CEp)*C;

alfa=10000;

function [LME, LMI, OBJ]=AWBT_EBRL(XLIST)
    
[X,Q2,Th,Kh,Yh,gama_w,gama_d,gama_v]= XLIST(:)
LME=list(X-X',Q2-Q2')
LMI=list(-([2*(B*Th+Th'*B'+A*X+X'*A'),zeros(n,n),2*G,2*E,zeros(n,ny),X*C';zeros(n,n),alfa*eye(n,n)+A'*P1'*Q2'+Q2*P1*A+Yh*P2*A+A'*P2'*Yh'-Kh*C-C'*Kh',Q2*P1*G-Kh+Yh*P2*G,zeros(n,nd),Q2*H1+Yh*H2,zeros(n,ny);2*G',G'*P1'*Q2'-Kh'+G'*P2'*Yh',-gama_w*eye(ny,ny),zeros(ny,nd+ny),eye(ny,ny);2*E',zeros(nd,n+ny),-gama_d*eye(nd,nd),zeros(nd,ny+ny);zeros(ny,n),H1'*Q2'+H2'*Yh',zeros(ny,ny+nd),-gama_v*eye(ny,ny),zeros(ny,ny);C*X',zeros(ny,n),eye(ny,ny),zeros(ny,nd+ny),-eye(ny,ny)]),gama_w,gama_d,gama_v,X,Q2)
OBJ_V=Q2*P1*E+Yh*P2*E;
OBJ=[norm(OBJ_V,'inf')]
//OBJ=[]

endfunction

X0=eye(n,n);
Q20=eye(n,n);
Th0=rand(nu,n);
Kh0=rand(n,ny);
Yh0=rand(n,ny);
gama_w0=1;
gama_d0=1;
gama_v0=1;

Init_guess=list(X0,Q20,Th0,Kh0,Yh0,gama_w0,gama_d0,gama_v0);

Mbound=1e5;
abstol=1e-15;
nu=2;
maxiters=100;
reltol=1e-10;
options=[Mbound,abstol,nu,maxiters,reltol];

Ans_LMI=lmisolver(Init_guess,AWBT_EBRL,options);
//Ans_LMI=lmisolver(Init_guess,AWBT_EBRL);

X=Ans_LMI(1);
Q2=Ans_LMI(2);
Th=Ans_LMI(3);
Kh=Ans_LMI(4);
Yh=Ans_LMI(5);
gama_w=Ans_LMI(6);
gama_d=Ans_LMI(7);
gama_v=Ans_LMI(8);

Y=Q2\Yh;
P=P1+Y*P2;
H=H1+Y*H2;
J=(eye(n,n)+(-E*CEp+Y*(eye(ny,ny)-(C*E)*CEp))*C)*B;
K=Q2\Kh;
T=Th/X;
N=P*A-K*C;
L=K-N*H;

savematfile("scilab_results.mat", "P", "H","J", "K","T", "N","L", "-v7.3");
