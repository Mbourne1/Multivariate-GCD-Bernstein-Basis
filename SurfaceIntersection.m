function [] = SurfaceIntersection

fxy_Bern = ...
    [
    0   0   0   0   ;
    0   -5  -5  0   ;
    0   -5  -5  0   ;
    0   0   0   0
    ];
% gxy_Bern = ...
%     [
%     -3  -3  -3  -3  ;
%     -3  2  2  -3  ;
%     -3  2  2  -3  ;
%     -3  -3  -3  -3  ;
%     ];

gxy_Bern = ...
    [
    -2 -2  -2  -2  ;
     -2 -2  -2  -2  ;
     -2 -2  -2  -2  ;
     -2 -2  -2  -2  ;
    ];

fxy-gxy


plot_fxy_gxy(fxy_Bern,gxy_Bern)

o1(fxy_Bern,gxy_Bern,3,3)


end