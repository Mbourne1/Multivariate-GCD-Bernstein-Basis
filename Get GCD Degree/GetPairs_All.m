function [t1t2_Pairs] = GetPairs_All(m1,m2,n1,n2)
%
%
%


% The total of t1+t2 must be between t and 2t
t1t2_Pairs = [];
for t1 = min(m1,n1):-1:0;
    for t2 = min(m2,n2):-1:0;
        
        new_row = [t1 t2];
        t1t2_Pairs = [t1t2_Pairs ; new_row];
        
    end
end

end