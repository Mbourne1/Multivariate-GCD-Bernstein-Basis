function [] = CoefficientScaling(m1, m2, n1 , n2, k1, k2)
% Set all coefficients of f(x,y) equal to one, determine scaling effect due
% to binomial terms in the subresultant matrices
%
% CoefficientScaling(m1, m2, n1, n2, k1, k2)
%
% % Inputs
%
% m1 : (Int) Degree of the polynomail f(x,y) with respect to x
%
% m2 : (Int) Degree of the polynomial f(x,y) with respect to y
%
% n1 : (Int) Degree of the polynomail g(x,y) with respect to x
%
% n2 : (Int) Degree of the polynomail g(x,y) with respect to y
%
% k1 : (Int) Index of subresultant matrix 
%
% k2 : (Int) Index of subresultant matrix
%
% >> CoefficientScaling(5, 5, 11, 11, 3, 2) 



% Get vector of coefficients of f(x,y)
fxy_vec = ones(m1 + 1 * m2 + 1);
fxy = GetAsMatrix(fxy_vec, m1, m2);

% Initialise arrays
arrT = cell(m1 + 1, m2 + 1);
arrDT = cell(m1 + 1, m2 + 1);
arrDTQ = cell(m1 + 1, m2 + 1);
arrTQ = cell(m1 + 1, m2 + 1);
arrDTQ_DenominatorsRemoved = cell(m1 + 1, m2 + 1);
% Scaling Effect of T

% For each coefficient a_{i1, i2}
for i1 = 0 : 1 : m1
    for i2 = 0 : 1 : m2
        
        % For each column of the subresultant matrix partitions
        for j1 = 0 : 1 : n1 - k1
            for j2 = 0 : 1 : n2 - k2
                
                
                arrT{i1 + 1, i2 + 1}(j1 + 1, j2 + 1) = ...
                    nchoosek(m1, i1) ...
                    * ...
                    nchoosek(m2, i2);
                
                arrDT{i1 + 1, i2 + 1}(j1 + 1, j2 + 1) = ...
                    (...
                    nchoosek(m1, i1) ...
                    * ...
                    nchoosek(m2, i2) ...
                    )...
                    ./ ...
                    ( ...
                    nchoosek(m1 + n1 - k1, i1 + j1) ...
                    * ...
                    nchoosek(m2 + n2 - k2, i1 + j1) ....
                    );
                
                arrTQ{i1 + 1, i2 + 1}(j1 + 1, j2 + 1) = ...
                    nchoosek(m1, i1) ...
                    * ...
                    nchoosek(m2, i2) ...
                    * ...
                    nchoosek(n1 - k1, j1) ...
                    * ...
                    nchoosek(n2 - k2, j2) ;
                                
                arrDTQ{i1 + 1, i2 + 1}(j1 + 1, j2 + 1) = ...
                    (...
                    nchoosek(m1, i1) ...
                    * ...
                    nchoosek(m2, i2) ...
                    * ...
                    nchoosek(n1 - k1, j1) ...
                    * ...
                    nchoosek(n2 - k2, j2) ...
                    )...
                    ./ ...
                    ( ...
                    nchoosek(m1 + n1 - k1, i1 + j1) ...
                    * ...
                    nchoosek(m2 + n2 - k2, i2 + j2) ....
                    );
                
                arrDTQ_DenominatorsRemoved{i1 + 1, i2 + 1}(j1 + 1, j2 + 1) =...
                    nchoosek(i1 + j1, i1) ...
                    * ...
                    nchoosek(i2 + j2, i2) ...
                    * ...
                    nchoosek(m1 + n1 - k1 - i1 - j1 , m1 - i1 ) ...
                    * ...
                    nchoosek(m2 + n2 - k2 - i2 - j2 , m2 - i2 ) ;
                
            end
        end
    end  
end






PlotValues(arrT);
PlotValues(arrDT);
PlotValues(arrTQ);
PlotValues(arrDTQ);
PlotValues(arrDTQ_DenominatorsRemoved);


end





function [] = PlotValues(arrT)
%
% % Inputs
%
%
% arrT : (Array of Matrices) Array of subresultant matrices T_{k1, k2}

figure()


[nArrayRows, nArrayCols] = size(arrT);

m1 = nArrayRows - 1;
m2 = nArrayCols - 1;


hold on

% For each coefficient a_{i1, i2}
for i1 = 0 : 1 : m1
    for i2 = 0 : 1 : m2
        
        temp_vec = GetAsVector_Version1(arrT{i1 + 1, i2 + 1}, m1, m2);
        plot(log10(temp_vec));
        
    end
end

hold off


end


