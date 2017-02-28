function [t1t2_Pairs] = GetPairs_All_3Polys(m1, m2, n1, n2, o1, o2)
%
% % Inputs
%
% [m1, m2] : Degree of polynomial f(x,y) with respect to x and y
%
% [n1, n2] : Degree of polynomial g(x,y) with respect to x and y
%
% [o1, o2] : Degree of polynomial h(x,y) with respect to x and y
%


% The total of t1+t2 must be between t and 2t
t1t2_Pairs = [];

for t1 = min([m1,n1,o1]):-1:0
    for t2 = min([m2,n2,o2]):-1:0
        
        new_row = [t1 t2];
        
        t1t2_Pairs = [t1t2_Pairs ; new_row];
        
    end
end

end