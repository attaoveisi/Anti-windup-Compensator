clc
clear
load('freq_vec.mat')
load('iddata_act1.mat')
load('iddata_act2.mat')
load('iddata_shaker.mat')

iddata_act1 = iddata_act1';
iddata_act2 = iddata_act2';
iddata_shaker = iddata_shaker';

ny = 1;
nu = 3;
nf = 8193;

iddata = zeros(ny,nu,nf);

for i = 1:nf
    iddata(:,1,i) = iddata_act1(:,i);
    iddata(:,2,i) = iddata_act2(:,i);
    iddata(:,3,i) = iddata_shaker(:,i);
end
    
freq_vec = freq_vec*2*pi;
