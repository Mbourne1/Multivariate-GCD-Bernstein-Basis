function [] = CoefficientScaling_WithCoefficients(m1, m2, n1 , n2, k1, k2)
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


fxy_vec = ones(m1 + 1 * m2 + 1);
fxy = GetAsMatrix(fxy_vec, m1, m2);

arrT = cell(m1 + 1, m2 + 1);
arrDT = cell(m1 + 1, m2 + 1);
arrDTQ = cell(m1 + 1, m2 + 1);
arrTQ = cell(m1 + 1, m2 + 1);
arrDTQ_DenominatorsRemoved = cell(m1 + 1, m2 + 1);
% Scaling Effect of T

[fxy_exact, gxy_exact, dxy_exact, uxy_exact, vxy_exact, m, n, t_exact] = ...
    Examples_GCD_Bivariate_2Polys(ex_num);


for i1 = 0 : 1 : m1
    
    for i2 = 0 : 1 : m2

        for j1 = 0 : 1 : n1 - k1
            
            for j2 = 0 : 1 : n2 - k2
                
                
                arrT{i1+1, i2+1}(j1 + 1, j2 + 1) = ...
                    nchoosek(m1, i1) ...
                    * ...
                    nchoosek(m2, i2);
                
                arrDT{i1+1, i2+1}(j1 + 1, j2 + 1) = ...
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
                
                arrTQ{i1+1, i2+1}(j1 + 1, j2 + 1) = ...
                    nchoosek(m1, i1) ...
                    * ...
                    nchoosek(m2, i2) ...
                    * ...
                    nchoosek(n1 - k1, j1) ...
                    * ...
                    nchoosek(n2 - k2, j2) ;
                
                
                
                
                arrDTQ{i1+1, i2+1}(j1 + 1, j2 + 1) = ...
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
                
                arrDTQ_DenominatorsRemoved{i1+1, i2+1}(j1 + 1, j2 + 1) =...
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

figure()


[nArrayRows, nArrayCols] = size(arrT);

m1 = nArrayRows - 1;
m2 = nArrayCols - 1;


hold on
for i1 = 0:1:m1
    
    for i2 = 0 : 1 : m2
        
        
        
        temp_vec = GetAsVector_Version1(arrT{i1+1, i2+1}, m1, m2);
        plot(log10(temp_vec));
        
    end
end

hold off


end


