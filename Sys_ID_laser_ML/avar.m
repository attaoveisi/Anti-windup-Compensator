%  Monte Carlo simulation to look at asymptotic variability of least squares
%  estimates vs theory for  a range of model structures.  
%
%  This script can use either the MATLAB system identification toolbox or 
%  the University of Newcastle Identification Toolbox (UNIT)
%
%  By default it uses UNIT, which may be dowloaded at 
%
%  http://sigpromu.org/idtoolbox
%
%  but if the user prefers to use the MATLAB toolbox, simply change the variable 
%  UNIT below to equal zero

UNIT = 0;   % Default is UNIT, set UNIT=0 to use MATLAB toolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters of simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 1000;                    %  Length of data record in seconds  
monte_number = 100;          %  Number of monte carlo simulations
fs = 1.0;                    %  Sampling frequency
Ts = 1/fs;
Tdelay = 0;                  %  Time delay in true system 
Ndelay = floor(fs*Tdelay);
throwaway = floor(T*fs/10);  %  Discard some samples to get rid of ic's
N = floor(T*fs); 
var = 0.01;                  %  Measurement noise variance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  System the user wants to estimate - select by setting the variable "model" below
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system = 9; 
switch system
 case{0}  % Ninness system
  den = [10,11,1]; num = 1;
  [numd,dend] = c2dm(num,den,Ts,'zoh');
  xi = real(roots(dend)); % z-domain poles   
  dend = poly([0.95]); numd = sum(dend);
 case{1}  % First example of OE example of paper
  xi = [0.8*ones(1,3),0.7*ones(1,2),0.4*ones(1,2)]; % z-domain poles
  den = real(poly(log(xi)));
  num = real(poly([0.3*ones(1,3)])); num = num*den(end)/num(end);  
   [numd,dend] = c2dm(num,den,Ts,'zoh');
 case{2}  % Resonant system
  xi = [0.95*ones(1,3),0.75*exp(j*pi/3),0.75*exp(-j*pi/3),0.60*ones(1,2)]; % z-domain poles
  den =real(poly(log(xi))); 
  num = den(length(den));
  [numd,dend] = c2dm(num,den,Ts,'zoh');
 case{3}  % Ortho-basis system
  xi = 0.5*ones(1,8); % z- domain poles
  dend = real(poly(xi));
  numd = sqrt(0.75)*real(poly(1/0.5*ones(1,7)));
  numd = numd*sum(dend)/sum(numd);  
 case{4}  % Randomly chosen system
  xi = rand(1,8); % z- domain poles
  dend = real(poly(xi));
  numd = real(poly(2*(ones(1,7)-2*rand(1,7))));
  numd = numd*sum(dend)/sum(numd);  
 case{5}  % Mid order system
  xi = [0.99,0.6,0.7]; % z- domain poles
  dend = real(poly(xi));  
  numd = real(poly([0.8,0.9]));
  numd = numd*sum(dend)/sum(numd);  
 case{6}  % Low (First) order system
  xi = [0.95]; % z- domain poles
  dend = real(poly(xi));  
  numd=[0,sum(dend)]; 
 case{7}  % Mid Order Resonant system
  xi = [0.75*exp(j*pi/3),0.75*exp(-j*pi/3),0.95*exp(j*pi/12),0.95*exp(-j*pi/12)]; % z-domain poles
  den = real(poly(log(xi))); 
  num = den(length(den));
  [numd,dend] = c2dm(num,den,Ts,'zoh');
 case{8}  % Another Low Order Resonant system
  xi = [0.95*exp(j*pi/12),0.95*exp(-j*pi/12)]; % z-domain poles
  den =real(poly(log(xi))); 
  num = den(length(den));
  [numd,dend] = c2dm(num,den,Ts,'zoh');
 case{9}  % Åström system
  numd = [1,0.5]; dend = [1,-1.5,0.7];
  xi = roots(dend); % z-domain poles
end;

% System in object oriented form
Gsys = tf(numd,dend,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Below here should not need any changing by the user
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Specify options about how estimation algorithms should run
OPT.n=throwaway; OPT.dsp=1;  OPT.fast = 1;

% Specify OE model structure and x-axis for frequency response estimation 
clear M;  w = logspace(log10(2*pi/(1000*Ts)),log10(pi/Ts),500);

% UNIT model structure initialisation
M.B=numd; M.A=dend; M.C=[]; M.D=[]; 
M.w = w; M.T = Ts;   

% Matlab ID toolbox model structure initialisation
m = idpoly(1,numd,1,1,dend,var,Ts);

%  Calculate the true frequency response.
ww = exp(j*w*Ts);
Gtrue = polyval(numd,ww)./polyval(dend,ww);

% Don't forget to include time delays
Gtrue = Gtrue .* exp(-j*w*Tdelay);

%  Now do Monte-Carlo over possible data realisations.
Gmonte = zeros(length(ww),monte_number);

for k=1:monte_number
k
  % Generate input excitation for reference
  r = randn(1,N); 

  % Simulate OE system
  noise = sqrt(var)*randn(size(r));
  u = r; y = filter(numd,dend,r) + noise; 

  % Estimate OE structure
  
  if UNIT
   Z.y=y; Z.u=u; 
   Gest = est(Z,M,OPT);     
   G = Gest.G(:);     
  else
   z=iddata(y(:),u(:),Ts);
   th = oe(z,m,'Trace','On');
   [x1,Bhat,x2,x3,Ahat] = th2poly(th);  
   G = polyval(Bhat,ww)./polyval(Ahat,ww);  
  end;
  
  Gmonte(:,k) = G;  
  
end; % End of loop over iterations

%  Calculate sample means and variances
for k=1:length(ww)
 meanG(k) = mean(Gmonte(k,:));
 varG(k)  = cov(Gmonte(k,:));
end;
  
%  Calculate Kn(w) = \sum_{k=0}^{p-1}|B_k|^2
xi = roots(dend);
% Now compute Kn(w) using these zeros.
Kn = zeros(length(w),1); 
for k=1:length(xi);
 Kn = Kn + (1-abs(xi(k))^2)./(abs(ww(:)-xi(k)).^2);
end;

% Compute variance approximation on input output dynamics

new = 1/(N-throwaway)*2*Kn(:)*var;  %Brett and Hakan's expression
old = 1/(N-throwaway)*(length(dend)-1)*var*ones(size(Kn)); % Lennart's expression

% Display results 

figure(1)
h1 = semilogx(w,10*log10(abs(varG(:))),'-g',w,10*log10(new),'--b',w,10*log10(old),'-.m');
legend('Sample Variability','New Theory','Exisiting Theory');
set(h1,'Linewidth',2); 
grid;
title('Variance of OE estimate of G vs existing and extended Theory')
xlabel('Frequency (normalised)')
ylabel('Variance (dB)')




