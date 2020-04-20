%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Script: atta_sim.m
%
% Author: MC Turner
%         Dept. of Engineering
%         Uni. of Leicester
%
% Date: 10th February 2006
%
% Modified: Adapted to model of Atta Oveisi
%           
%
% Purpose: To build Atta's plant, call anti-windup design functions and to
% run simulation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct plant

load identified_t
B = B/1e4;

Ap  =  A;
Bp  =  B;
Bpd =  H
Cp  =  C;
Dp  = zeros(1,2);
Dpd = zeros(1,1);


% Construct controller

%%%%%%% Kalma Filter Design
sysk=ss(Ap,[Bp Bpd],Cp,[Dp Dpd]);

[kest,L,P] = kalman(sysk,1e3,1,0);

%%%%%%% Optimal Gain
syst=ss(A,B,C,0);
[K,S,e] = lqry(syst,1e0,1e0*eye(2),0);

%%%%%%% Controller
rlqg = lqgreg(kest,K);

Ac = rlqg.A;
Bc = rlqg.B;
Cc = rlqg.C;
Dc = rlqg.D;


% Control limit

ubar = 0.5;

%%
%--------------------------------------------------------------------------

%--------------------------------
%
% Design full-order AW compensator
%
%--------------------------------

disp(' ');
disp('....Designing full-order AW compensator....');
disp(' ');

G   =   ss(Ap,Bp,Cp,Dp);% Plant without disturbance inputs
Wp  =   1e-1;
Wr  =   1*eye(1); % increasing robustness weight tends to give smaller poles

                           
[AWfull,gam] = fullorder_ctf(G,Wp,Wr); 

%%
%--------------------------------------------------------------------------
% Setup simulations

tstop = 14;    % Stop time
tstep = 1e-4  % Time-step for solver

[np,m]  = size(Bp);
[np,nd] = size(Bpd);
[nc,p]  = size(Bc);


% Simulate linear model
baru   = inf;
AW     = ss([],zeros(0,m),zeros(m+p,0),zeros(m+p,m));
     
sim('genericAW_atta');
      
figure(1);
subplot(311);
plot(y.time,y.signals.values,'b');
ylabel('Output response');
title('Simulation');
hold on;
subplot(312);
ul=plot(um.time,um.signals.values(:,1),'b');
ylabel('Control response (1)');
hold on;
subplot(313);
ul=plot(um.time,um.signals.values(:,2),'b');
hold on;
ylabel('Control response (2)');
xlabel('Time [sec]');
hold on;


% Simulate saturated model -no AW

baru   = ubar;
AW     = ss([],zeros(0,m),zeros(m+p,0),zeros(m+p,m));
     
sim('genericAW_atta');
      
figure(1);
subplot(311);
plot(y.time,y.signals.values,'r');
subplot(312);
us=plot(um.time,um.signals.values(:,1),'r');
subplot(313);
us=plot(um.time,um.signals.values(:,2),'r');


% Simulate saturated model - with full-order AW

baru   = ubar;
AW     = AWfull;
     
sim('genericAW_atta');
      
figure(1);
subplot(311);
plot(y.time,y.signals.values,'g');
axis([0 tstop -12 12]); 
subplot(312);
uaw=plot(um.time,um.signals.values(:,1),'g');
subplot(313);
uaw=plot(um.time,um.signals.values(:,2),'g');
%axis([0 tstop -60 60]); 

%legend([ul(1),us(1),uaw(1)],'Linear','Sat w/o AW','Sat w/AW');

% hline = findobj(gcf, 'type', 'line');
% set(hline,'LineWidth',2);

