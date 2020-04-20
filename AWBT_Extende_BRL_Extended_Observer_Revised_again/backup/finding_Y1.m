function output = finding_Y(Y)

E = [0;-5.42650000000000;0;5.07950000000000;0;5.00350000000000];
C = [0,-5.10120000000000,0,4.02600000000000,0,3.27260000000000];

% E=E(2:5,:);
% C=C(:,2:5);

dum3 = size(C);
ny = dum3(1,1);
CEp=(((C*E)'*(C*E))^-1)*(C*E)';
P1 = eye(6,6)-E*CEp*C;

% P1 = eye(4,4)-E*CEp*C;

P2 = (eye(ny,ny)-(C*E)*CEp)*C;

if ny == 2
    Y1 = Y(1,1:6)';
    Y2 = Y(1,7:end)';
%     Y1 = Y(1,1:4)';
%     Y2 = Y(1,5:end)';
    Y = [Y1 Y2];
else
    Y = Y';
end

output1 = (P1+Y*P2)*E;
output = sumsqr(output1);
