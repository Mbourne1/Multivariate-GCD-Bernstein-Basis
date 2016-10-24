
function [f] = GetSymbolicPoly(factors_f)
% Given the array of factors of f(x,y) 
%
% Inputs
%
% factors_f : Factors of polynomial f(x,y)

syms x y

f = factors_f{1};
for i = 2:1:length(factors_f)
   f = f * factors_f{i}; 
end



end