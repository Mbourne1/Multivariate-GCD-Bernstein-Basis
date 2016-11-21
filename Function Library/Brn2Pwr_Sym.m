function [coef_pwr_sym] = Brn2Pwr_Sym(Coeff)

% Given the series of control points generate the power polynomial

[~,c] = size(Coeff);
deg_f = c-1;

t = sym('t');

sum =0;

for i =0:1:deg_f
    expand(Coeff(i+1) .* nchoosek(deg_f,i) .* ((1-t).^(deg_f-i)) .* (t.^i));
    sum = sum + Coeff(i+1) .* nchoosek(deg_f,i) .* ((1-t).^(deg_f-i)) .* (t.^i);
    
end

coef_pwr_sym = expand(sum);

end