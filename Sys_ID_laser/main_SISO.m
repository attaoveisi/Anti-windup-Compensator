clc
clear
close all

%% Fast Fourrier Transformation parameters
N = 6400; % 50,100,200,400,800,1600,3200,6400
nFFT = 2.56*N;

% Baseband
fmax = 800;  % frequency span (Hz)
fs = 2.56*fmax; % Nyquist frequency (Hz)

dt = 1/fs;  % sampling time (sec)
df = fmax/nFFT;  % sampling frequency (Hz)
frequency_vector = df:df:fmax;
% frequency_vector = fmax*linspace(0,1,nFFT/2+1);

T_total = nFFT*dt; % duration of the experiment (sec)
time_vector = 0:dt:T_total-dt;

%% Input signal generation with window and overlap
n_inp = 1; % number of input channels
n_out = 1; % number of output channels

% Windowing (Hanning window): w(n)=0.5(1?cos(2?n/N)), 0?n?N L=N+1
% w = hann(L,'sflag') returns an L-point Hann window using the window sampling specified by 'sflag',
% which can be either 'periodic' or 'symmetric' (the default). 
% The 'periodic' flag is useful for DFT/FFT purposes, such as in spectral analysis.

hann_win_inp = zeros(nFFT,n_inp);
for i = 1:n_inp
    hann_win_inp(:,i) = hann(nFFT,'periodic');
end

hann_win_out = zeros(nFFT,n_inp);
for i = 1:n_out
    hann_win_out(:,i) = hann(nFFT,'periodic');
end

% averaging and overlap properties
n_aver = 100; % number of averaging

% overlap parameters
overlap = 75; 
r_overlap = overlap/100; % rate of overlap (%)

T_total = n_aver*T_total-(n_aver-1)*r_overlap*T_total; % total duration of input signal after windowing and averaging (sec)
total_samples = round(n_aver*nFFT-(n_aver-1)*r_overlap*nFFT); % total number of samples
total_time_vector = 0:dt:T_total-dt;
[a1,a2] = size(total_time_vector);
% Creating the input signal (uniformly distributed random signal)
% load act1
% data3 = act1;
% 
% inp = data3.Y(1).Data;

load act2
data3 = act2;

inp = data3.Y(2).Data;

% load shaker
% data3 = shaker;
% 
% inp = data3.Y(4).Data;

inp = inp(1,end-(a2-1):end)';

% load act1_only
% inp = data2.Y(1).Data;

% figure(1)
% subplot(2,1,1)
% plot(total_time_vector,inp)
% xlabel('time (sec)')
% ylabel('amplitude')
% subplot(2,1,1)
% n_bins = 500; % number of Histrogram bins
% subplot(2,1,2)
% title('Simulated Histogram of White uniform Noise');
% for i = 1:n_inp
%     [PDF_inp,bins] = hist(inp(:,i),n_bins);
%     hold on
%     bar(bins,PDF_inp/trapz(bins,PDF_inp));
% end
% hold off
% xlabel('bins')
% ylabel('PDF of input signals')
% plot(total_time_vector,inp)

%% Output analysis
out = data3.Y(3).Data;
% out = data3.Y(4).Data;

% out = data2.Y(3).Data;
% out = data2.Y(4).Data;


out = out(:,end-(a2-1):end)';

if overlap > 0
    n_samples_in_Overlap = floor((overlap*nFFT) / 100); 
    nFrames = n_aver; 
else
    n_samples_in_Overlap= nFFT;
    nFrames=floor(total_samples/n_samples_in_Overlap)-1;
end

S_inp_out = zeros(nFFT/2+1,n_inp*n_out);
k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        [S_inp_out_dum, freq_after_overlap] = cpsd(inp(:,i),out(:,j),hann_win_inp(:,1),overlap,nFFT,df);
        S_inp_out(:,k) = S_inp_out_dum;
    end
end

S_inp_inp = zeros(nFFT/2+1,n_inp);
for i = 1:n_inp
    [S_inp_inp(:,i), ~] = pwelch(inp(:,i),hann_win_inp(:,i),overlap,nFFT,df);
end

S_out_out = zeros(nFFT/2+1,n_out);
for j = 1:n_out
    [S_out_out(:,j), ~] = pwelch(out(:,j),hann_win_out(:,j),overlap,nFFT,df);
end

k = 0;
for i = 1:n_inp
    for j = 1:n_out
        k = k+1;
        H_2(:,k) = S_out_out(:,j)./S_inp_out(:,k);
    end
end

% figure(3)
% k = 0;
% for i = 1:n_inp
%     for j = 1:n_out
%         k = k+1;
%         subplot(n_inp,n_out,k)
%         semilogy(freq_after_overlap*1/dt/df,abs(H_2(:,k)))
%         xlim([0 800])
%     end
% end

%% SISO ID
ww=(2*pi*freq_after_overlap*1/dt/df)'; ff=ww/(2*pi);
j=sqrt(-1); 
B_Complex=H_2;  %if acceleratoin is measured (FRF system az inja shuru mikonam man)
Gjw=B_Complex;   %if displacement is measured (chon acc andaze giri shode vali man lazem nadaram)
% =========================================================================
% transform continuous data to discrete using bilinear tansfrmation
% !!!!!!!!!!!!!!!! see help bilinear
select_w=1; % chon ferekans payin javabesh bade arc tan dar nazdike sefr javabesh khub nis.
T=5/max(ww);
wwd=2*atan(T*ww(select_w:end)/2);
wwd=wwd';
wwd;
Gzks=Gjw(select_w:end);
% Gzks(1,1)=500+500j;      % because Gjw(1,1) is inf         
% =========================================================================
m=1;                 % number of input
p=1;                 % number of output   
[M,MM]=size(Gzks);    % number of data
j=sqrt(-1);
n=16;                % order of system
q=n*5;        % q should be greater than system order (Eq 6)


% ===================== step 1: compute matrix G  =========================
% **************  according to equation (47) or (10) in conf paper*********
% G
for i=1:q
    for k=1:M
    Y_qp_mM(p*(i-1)+1:p*i,m*(k-1)+1:m*k)=exp(j*(i-1)*wwd(k,1))*Gzks((k-1)*p+1:k*p,:);
    end
end

% ==================== step 1: compute matrix Wm  =========================
% *********************  according to equation (48) ***********************
% Wm
for i=1:q
    for k=1:M
    U_m(p*(i-1)+1:p*i,m*(k-1)+1:m*k)=exp(j*(i-1)*wwd(k,1))*eye(m,m);
    end
end
U_re=[real(U_m),imag(U_m)];
Y_re=[real(Y_qp_mM),imag(Y_qp_mM)];
% ================= step 2: compute QR factorizatin   =====================
% ******* according to equation (62) using lqdecomposition function *******
[Q,R]=qr([U_re' Y_re'],0);
R22=R(end-q+1:end,end-q+1:end);
% ==================== step 3: compute SVD  ===============================
% ******* according to equation (63) using lqdecomposition function *******
[U_hat,S_hat,V_hat] = svd(R22');

% =================== step 4: detemine system order  ======================
% **** according to equation (64) and estimate of observability matrix ****
U_hat_s=U_hat(:,1:n);

% ================= step 5: detemine A_hat & C_hat  =======================
% ***************** according to equations (65) & (66) ********************
O_hat_upperline=U_hat_s(p+1:end,:);
O_hat_underline=U_hat_s(1:end-p,:);
A_hat=inv(O_hat_underline'*O_hat_underline)*O_hat_underline'*O_hat_upperline;
C_hat=U_hat_s(1:p,:);
% ================= step 6: detemine B_hat & D_hat  =======================
% ********************* according to equation (67) ************************
teta(1:n+1,1)=0;       % the vector contains of B_hat and D_hat coefficient Az least square problem
for i=1:(M)
    C_hattimeinvA_hat(i,:)=C_hat*inv((exp(j*wwd(i,1)))*eye(size(A_hat))-A_hat);
end
zaribD(1:M,1)=1;
B_hat=teta(1:n,1);
D_hat=teta(n+1,1);
epsilon=Gzks(1:M)-D_hat*C_hattimeinvA_hat*B_hat; % estimation error
norm(epsilon)
jakoobiyan=[C_hattimeinvA_hat zaribD];
jakoobiyan_re=[real(jakoobiyan);imag(jakoobiyan)];
% teta=inv((jakoobiyan_re'*jakoobiyan_re))*((jakoobiyan_re')*[real(Gzks(1:M));imag(Gzks(1:M))]);
teta=inv((jakoobiyan_re'*jakoobiyan_re))*((jakoobiyan_re')*[real(Gzks(1:M));imag(Gzks(1:M))]);
B_hat=teta(1:n,1);
D_hat=teta(n+1,1);
epsilon=Gzks(1:M)-D_hat*C_hattimeinvA_hat*B_hat;
norm(epsilon)
% =============== step 7: estimated transfer function  ====================
% ********************* according to equation (68) ************************
G_hat=ss(A_hat,B_hat,C_hat,D_hat,T);
% figure
% bode(G_hat)

% % % % % % % % % % % % % % % % % % % % % %   
sysc=d2c(G_hat,'tustin');
figure
plot(ff(select_w:end),20*log10(abs(Gjw(select_w:end))));
hold on
[magG,phaseG,w]=bode(sysc,10.5*2*pi:pi:2*pi*800);
  for k=1:size(w)
      bb(k,1)=magG(:,:,k);
  end
  plot(w/2/pi,20*log10(bb),'r')
  
%   xlim([0 800])


