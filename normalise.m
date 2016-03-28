function [fxy_n] = normalise(fxy)
% Normalise coefficients of the matrix by dividing by the first entry.

fxy_n = fxy./fxy(1,1);


end
