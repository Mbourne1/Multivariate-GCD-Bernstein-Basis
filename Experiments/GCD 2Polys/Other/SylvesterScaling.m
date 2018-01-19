function [] = SylvesterScaling(m1, m2, n1, n2, k1, k2) 
%
% % Inputs
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
%
% n1 : (Int) Degree of g(x,y) with respect to x
%
% n2 : (Int) Degree of g(x,y) with respect to y
%
% k1 : (Int) Degree of d(x,y) with respect to x
%
% k2 : (Int) Degree of d(x,y) with respect to y
%
%
% % Examples 
%
% >> SylvesterScaling(10, 15, 9, 6, 2, 1)
%
% >> SylvesterScaling(5, 4, 7, 6, 2, 2)


close all; clc;


arrMethods = {'T', 'DT', 'TQ', 'DTQ', 'DTQ Denominator Removed'};

for i = 1 : 1 : length(arrMethods)
    
    method_name = arrMethods{i};
   
    Sk = BuildSylvesterMatrix(m1, m2, n1, n2, k1, k2, method_name);
    
    PlotHeatMap(Sk, method_name)
    
    
end



end



function Sk = BuildSylvesterMatrix(m1, m2, n1, n2, k1, k2, method_name)
%
% % Inputs
%
% m1 : (Int)
%
% m2 : (Int)
%
% n1 : (Int)
%
% n2 : (Int)
%
% k1 : (Int)
%
% k2 : (Int)
%
% method_name : (String)


    fxy = ones(m1 + 1, m2 + 1);
    gxy = ones(n1 + 1, n2 + 1);

    global SETTINGS
    SETTINGS.SYLVESTER_MATRIX_VARIANT = method_name;
    SETTINGS.VECTORISATION_METHOD = 'Version 1';
    Sk = BuildSubresultant_Bivariate_2Polys(fxy, gxy, k1, k2);

end

function [] = PlotHeatMap(Sk, method_name)
%
% % Inputs
%
% Sk : (Matrix)
%
% method_name : (String)

Sk = abs(Sk);
Sk = log(Sk);


figure('name',method_name)

h = heatmap(Sk);
%set(h,'GridVisible','off')
    
    
    



end