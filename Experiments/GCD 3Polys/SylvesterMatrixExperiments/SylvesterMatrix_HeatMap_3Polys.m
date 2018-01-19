function [] = SylvesterMatrix_HeatMap_3Polys(m1, m2, n1, n2, o1, o2, k1, k2)
% SylvesterMatrix_HeatMap_3Polys(m1, m2, n1, n2, o1, o2, k1, k2)
%
%
% m1 : (Int) Degree of f(x,y) with respect to x
%
% m2 : (Int) Degree of f(x,y) with respect to y
%
% n1 : (Int) Degree of g(x,y) with respect to x
%
% n2 : (Int) Degree of g(x,y) with respect to y
%
% o1 : (Int) Degree of h(x,y) with respect to x
%
% o2 : (Int) Degree of h(x,y) with respect to y
%
% k1 : (Int) Degree of d(x,y) with respect to x
%
% k2 : (Int) Degree of d(x,y) with respect to y
%
%
% >> SylvesterMatrix_HeatMap_3Polys(5, 5, 4, 3, 7, 9, 1, 2)

% Initialise polynomails f(x,y), g(x,y) and h(x,y)
fxy = ones(m1 + 1, m2 + 1);
gxy = ones(n1 + 1, n2 + 1);
hxy = ones(o1 + 1, o2 + 1);


% Get an array of the subresultant matrix variants
arrSylvesterFormats = {'T', 'DT', 'TQ', 'DTQ'};

global SETTINGS
SETTINGS.VECTORISATION_METHOD = 'Version 1';
SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS = '2';

% For each of the subresultant matrix variants
for i = 1 : 1 : length(arrSylvesterFormats)

    
    SETTINGS.SYLVESTER_MATRIX_VARIANT = arrSylvesterFormats{i};
    
    % Build the k1,k2-th subresultant matrix
    Sk = BuildSubresultant_Bivariate_3Polys(fxy, gxy, hxy, k1, k2);

    figure_name = sprintf('%s : Heat Map', SETTINGS.SYLVESTER_MATRIX_VARIANT);
    figure('Name',figure_name)
    
    Sk_rounded = log10(abs(Sk));
    
    hm = heatmap((Sk_rounded));
    %colormap(hm, 'hot')
    
    
end




end