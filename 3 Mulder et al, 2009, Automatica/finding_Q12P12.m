function output = finding_Q12P12(Y_dum)

[~,n] = size(Y_dum);
m = sqrt(n/2);
np = m;
global D_hat B_hat P11 Bp C_hat Cy Q11 A_hat Ap

Q12 = zeros(m,m);
P12 = zeros(m,m);
k = 1;
for i = 1:m
    for j = 1:m
        Q12(i,j) = Y_dum(1,k);
        k = k+1;
    end
end
for i = 1:m
    for j = 1:m
        P12(i,j) = Y_dum(1,k);
        k = k+1;
    end
end
% Q12 = rand(size(Q12));
% P12 = rand(size(P12));
output1 = max(max(abs(Q12*P12'-(eye(np,np)-Q11*P11))));

if  0 %(cond(P12) > 1e9 || cond(Q12') > 1e9)
    output = 1e6;
else
    Dk = D_hat;
    Bk = inv(P12)*(B_hat-P11*Bp*Dk);
    Ck = (C_hat-Dk*Cy*Q11)*inv(Q12');
    Ak = inv(P12)*(A_hat-P12*Bk*Cy*Q11-P11*Bp*Ck*Q12'-P11*(Ap+Bp*Dk*Cy)*Q11)*inv(Q12');

    if ( sum(sum(Ak(:,:) == inf)) || sum(sum(isnan(Ak(:,:)))) )
        dum1 = 10;
    else
        dum1 = max(real(eig(Ak)));
    end

    if dum1 >= 0
        output2 = 10;
    else
        output2 = 0;
    end

    output = output1 + output2;
end

