function [] = SylvesterMatrix_HeatMap_3Polys(m1, m2, n1, n2, o1, o2, k1, k2)
% SylvesterMatrix_HeatMap_3Polys(m1, m2, n1, n2, o1, o2, k1, k2)
%
%
%
%
%
%
%
% >> SylvesterMatrix_HeatMap_3Polys(5, 5, 4, 3, 7, 9, 1, 2)

fxy = ones(m1 + 1, m2 + 1);
gxy = ones(n1 + 1, n2 + 1);
hxy = ones(o1 + 1, o2 + 1);



arrSylvesterFormats = {'T','DT','TQ','DTQ'};

global SETTINGS
SETTINGS.VECTORISATION_METHOD = 'Version 1';
SETTINGS.SYLVESTER_MATRIX_3POLY_N_EQUATIONS = '2';


for i = 1 : 1 : length(arrSylvesterFormats)

    SETTINGS.SYLVESTER_BUILD_METHOD = arrSylvesterFormats{i};
    
    Sk = BuildSubresultant_Bivariate_3Polys(fxy,gxy,hxy,k1,k2);

    
    figure_name = sprintf('%s : Heat Map', SETTINGS.SYLVESTER_BUILD_METHOD);
    figure('Name',figure_name)
    
    Sk_rounded = log10(abs(Sk));
    
    hm = heatmap((Sk_rounded));
    %colormap(hm, 'hot')
    
    
end




end