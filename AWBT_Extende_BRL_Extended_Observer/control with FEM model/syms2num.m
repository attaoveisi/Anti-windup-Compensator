function[TFnum] = syms2num(G)
[symNum,~] = numden(G); %Get num and den of Symbolic TF
TFnum = sym2poly(symNum);    %Convert Symbolic num to polynomial