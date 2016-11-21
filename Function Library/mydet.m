% NOT WORKING

function [d] = mydet(mat)

% Get the determinant of the matrix whose entries contain x y and float
% parts.

end


function [d] = mydet_sub(mat)
% Given a two by two matrix with three entries, 

% a b
% c d

% ad - bc

%(a(x) + a(y) + a(z)) * (d(x) + d(y)+d(z))

mymult(a,d)

end

function ans = mymult(a,d)

a = [0 a(1) a(2) a(3)];
d = [0 d(1) d(2) d(3)];

C = BuildC(a)

ans = C*d


end
