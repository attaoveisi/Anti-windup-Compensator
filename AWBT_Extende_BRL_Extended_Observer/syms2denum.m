function[TFden] = syms2denum(G)
[~,symDen] = numden(G); %Get num and den of Symbolic TF 
TFden = sym2poly(symDen);    %Convert Symbolic den to polynomial