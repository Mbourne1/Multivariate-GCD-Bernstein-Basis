function [t1t2_Pairs] = GetPairs(m,m1,m2,n,n1,n2)
%
%
%


method = 'All';

switch method
    case 'All'  % Use all possible (t1,t2) combinations
        % The total of t1+t2 must be between t and 2t
        t1t2_Pairs = [];
        for t1 = min(m1,n1):-1:0;
            for t2 = min(m2,n2):-1:0;
                
                new_row = [t1 t2];
                t1t2_Pairs = [t1t2_Pairs ; new_row];
                
            end
        end
        
    case 'Refined' % Use only a subset of possible (t1,t2) combinations.
        % The total of t1+t2 must be between t and 2t
        t1t2_Pairs = [];
        
        for t1 = t:-1:0;
            for t2 = t:-1:0;
                
                condA = n1-t1 + n2 -t2 >= n-t;
                condB = m1-t1 + m2 -t2 >= m-t;
                condC = t1 <= n1 && t1 <= m1;
                condD = t2 <= n2 && t2 <= m2;
                condE = n1 - t1 + n2 - t2 <= 2*(n-t);
                condF = m1 - t1 + m2 - t2 <= 2*(m-t);
                condG = n1 - t1 <= n-t;
                condH = n2 - t2 <= n-t;
                condI = m1 - t1 <= m-t;
                condJ = m2 - t2 <= m-t;
                
                
                if condA && condB && condC && condD && condE && condF...
                        && condG && condH && condI && condJ
                    new_row = [t1 t2];
                    t1t2_Pairs = [t1t2_Pairs ; new_row];
                end
                
            end
        end
end

end