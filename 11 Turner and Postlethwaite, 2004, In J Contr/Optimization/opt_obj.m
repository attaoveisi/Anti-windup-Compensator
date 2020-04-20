function output = opt_obj(vars)
global Ap Bp Hp Cp Dp Ak Bk Ck Dk Af Bf Hf Cf Df theta F1 F2
nu = 2;
ny = 2;
F1 = tf(vars(1)*2*pi,[1 vars(2)*2*pi])*eye(nu);
F2 = tf(vars(3)*2*pi,[1 vars(4)*2*pi])*eye(ny);

options = simset('SrcWorkspace','current');
sim('Opt_AW',[],options)
output = sum((output.signals.values(:,1)).^2)+sum((output.signals.values(:,2)).^2);