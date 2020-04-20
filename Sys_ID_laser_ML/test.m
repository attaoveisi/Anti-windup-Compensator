clc
clear
close all
global dsp; if isempty(dsp), clear global dsp; dsp=1; end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PZT sensor (y)
load('iddata.mat')

ii = 6400;

iddata_act1 = act1_1(1:ii,3)'+act1_1(1:ii,4)'*1i;
iddata_act2 = act2_1(1:ii,3)'+act2_1(1:ii,4)'*1i;
iddata_shaker = shaker_1(1:ii,3)'+shaker_1(1:ii,4)'*1i;
freq_vec = act1_1(1:ii,2)*2*pi;

ny = 1;
nu = 2;
nf = ii;
N = nf;

iddata = zeros(ny,nu,nf);

for i = 1:nf
    Y(:,1,i) = iddata_act1(:,i);
    Y(:,2,i) = iddata_act2(:,i);
%     Y(:,3,i) = iddata_shaker(:,i);
end

n = 12;
m = nu;
p = ny;

w = real(freq_vec)';
Q=1*eye(n); 
R=1e-1*eye(p);   %Initial guess at covariances
T = 1/(800*2.56);

%%
z.y = Y; z.w = w;          % Specify measurements and frequencies they were obtained at
z.T = T*0;       % Specify sample time
m.nx = n;                  % Specify model order
m.type='ss';               % This is the default as well
m.w = w(:);
m.op = 's';                % Specify a continuous time model

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedure runs
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% opt.alg = 'sid';           % Specify subspace method
% opt.lag=round((N-10)/2);

% opt.alg = 'gn';            % Specify Gauss--Newton Search          
% opt.par   = 'ddlc'; 
% opt.cost = 'det'; 
% opt.cost = 'trace';
% opt.op = 's';
% opt.dir = 'trust';
% opt.ngt   = 0;


% opt.alg = 'em';            % Specify Continuous time state space model via EM Algorithm
% opt.optit   = 100;
% opt.stoptol = 1e-5; 

Ms.w = w; 
Ms.nx = n; 
Ms.op = 's'; 
Ms.T = 0; 
Ms.type='ss';               % This is the default as well
opt.alg = 'em';            % Specify Gauss--Newton Search


%% Global parameters
% opt.dsp = dsp;
% opt.miter = 1000; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% identified_plant_ss = est(z,m,opt);     % Subspace-based estimate
identified_plant_ss = fsid(z,Ms,opt);

%%
gt = identified_plant_ss;
A = gt.ss.A;
B = gt.ss.B(:,1);
H  = gt.ss.B(:,2);
C  = gt.ss.C;

[n,n_inp] = size(B);
[~,nw] = size(H);
[n_out,~] = size(C);

D = zeros(n_out,n_inp+nw);

[gt_reduced,redinfo1] = balancmr(ss(A,[B H],C,D),10); % display Hankel SV plot

rank(obsv(A,C))
rank(ctrb(A,B))
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dsp 
     data.G = Y; 
     data.w = w(:); 
     data.disp.legend = 'Data';
     data.disp.unit = 'hz';
     data.disp.colour    = 'b';
     data.disp.linestyle = '-';
     data.disp.linewidth = 0.5;
     
     gt.disp.legend = 'identified';
     gt.disp.colour    = 'r';
     gt.disp.linestyle = '-';
     gt.disp.linewidth = 0.5;

     showbode(data,gt);

end

figure;
bode(ss(A,[B H],C,D))
hold on
bode(gt_reduced)

rank(ctrb(gt_reduced.A,gt_reduced.B))
rank(obsv(gt_reduced.A,gt_reduced.C))