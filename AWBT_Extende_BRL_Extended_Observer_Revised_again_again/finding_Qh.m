function output = finding_Qh(Qh1)

B = [-8.88712865173690,3.42023387272080;4.60565774420267,-2.14780344563457;0.685686233358913,0.890353966519012;5.79297575782585,3.95333481790732];

nu = 2;

global Q1

j = 1;
for i = 1:nu
    Qh(j,:) = Qh1((j-1)*nu+1:j*nu);
    j = j+1;
end

output1 = Q1*B-B*Qh;
output = sumabs(output1);
