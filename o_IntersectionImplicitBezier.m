function [] = o_IntersectionImplicitBezier(ex_num,el,bool_preproc, bool_sntln,seed)
% Get the points of intersection between the explicitly defined surface
% z = f(x,y) and the the Bézier surface patch

global BOOL_Q
BOOL_Q = 'y';

global BOOL_PREPROC
BOOL_PREPROC = bool_preproc;

global BOOL_NOISE
BOOL_NOISE = 'y';

global bool_degreeMethod
bool_degreeMethod = '2';

global bool_remove_cols
bool_remove_cols = 'n';

global BOOL_SNTLN
BOOL_SNTLN = bool_sntln;

global PLOT_GRAPHS
PLOT_GRAPHS = 'n';

global SEED
SEED = seed;

% Get the implicit surface
fxy = ...
    [
        0 1;
        1 0;
    ];

%% Get the parametric Bezier surface

% %    P(:,:,1): x-coordates of control points as 4 x 4 matrix 
% %    P(:,:,2): y-coordates of control points as 4 x 4 matrix 
% %    P(:,:,3): z-coordates of control points as 4 x 4 matrix

P(:,:,1) = [0 0 0 0; 2 2 2 2; 3 3 3 3; 4 4 4 4];
P(:,:,2) = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4];
P(:,:,3) = [1 5 7 5; 1 5 7 5; 1 4 3 4; 4 2 3 4];

Q = bezierpatchinterp(P);
%plotbezierpatch3D(P,Q);

az=21;  %azimuth
el=19;  %elevation. 
figure
hold on

[Y,X,Z] =  ndgrid(linspace(-10,10,10),linspace(-10,10,10),linspace(-10,10,10));
V = X - 2*Y + Z; % evaluate your implicit function over the N-D grid
p = patch(isosurface(X,Y,Z,V,0));
isonormals(X,Y,Z,V,p);
set(p,'FaceColor','b','EdgeColor','k','FaceAlpha',0.5) % 'EdgeColor','none'
daspect([1 1 1])
axis square;
grid on;
camlight 
view(-27,46);
lighting gouraud

surface(Q(:,:,1),Q(:,:,2),Q(:,:,3),'FaceColor','green')
title('\bf Bezier Patch using Interpolated Points');
view(3);
box;  
view(az,el)

% Get the implicit surface
x = sym('x');
y = sym('y');
z = sym('z');

C1 = x  + y - z;

% Substitute in values for the parametric
C2_x = P(:,:,1);
C2_y = P(:,:,2);
C2_z = P(:,:,3);

C1 = subs(C1,{x,y,z},{C2_x,C2_y,C2_z})

GetRoots(C1)




end