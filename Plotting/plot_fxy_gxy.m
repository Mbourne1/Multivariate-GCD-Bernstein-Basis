function [] = plot_fxy_gxy(fxy, gxy)
%
%
% % Inputs
% 
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% gxy : (Matrix) Coefficients of polynomial g(x,y)

%%  Get the dimensions of f(x,y)


% Plot the control points of the surface of f(x,y)
fxy_coordinates_CP = getCoordinates_CP(fxy);
gxy_coordinates_CP = getCoordinates_CP(gxy);

% Plot the control points of the surface of g(x,y)
% Plot the control points in a three dimensional scatter
figure('name','Bernstein Control Points', 'NumberTitle','off')
hold on
title('The Bernstein Control points of f(x,y) and g(x,y)')
scatter3(fxy_coordinates_CP(:,1), fxy_coordinates_CP(:,2), fxy_coordinates_CP(:,3),'Filled')
scatter3(gxy_coordinates_CP(:,1), gxy_coordinates_CP(:,2), gxy_coordinates_CP(:,3),'Filled')
hold off


% Plot 3
% % Plot the Bézier surface by evaluating at a series of points

% Get the surface coordinates
fxy_Surf_coord = get_Coordinates_Curve(fxy);
gxy_Surf_coord = get_Coordinates_Curve(gxy);


% % plot the exact Bernstein data and the noisy bernstein data
figure('name','Bernstein Surfaces','NumberTitle','off')
title('Bernstein Surface')
hold on
% Plot exact data
scatter3(fxy_Surf_coord(:,1),fxy_Surf_coord(:,2),fxy_Surf_coord(:,3),'Filled')
scatter3(gxy_Surf_coord(:,1),gxy_Surf_coord(:,2),gxy_Surf_coord(:,3),'Filled')
hold off

%% Plot 4
% % Plot the bezier surface from the points above

% Get x y and z components of the exact Bernstein data
fxy_x = fxy_Surf_coord(:,1);
fxy_y = fxy_Surf_coord(:,2);
fxy_z = fxy_Surf_coord(:,3);

gxy_x = gxy_Surf_coord(:,1);
gxy_y = gxy_Surf_coord(:,2);
gxy_z = gxy_Surf_coord(:,3);


% Get step size
steps = 0:.05:1;

% obtain a grid of points (x,y)
[f_XI,f_YI] = meshgrid(steps, steps);
[g_XI,g_YI] = meshgrid(steps, steps);



% now interpolate - find z values for the grid points
f_ZI = griddata(fxy_x,fxy_y,fxy_z,f_XI, f_YI);
g_ZI = griddata(gxy_x,gxy_y,gxy_z,g_XI, g_YI);

% % Plot the surface
figure(4)
title('Exact Bernstein Surface')
hold on
g1 = surf(f_XI,f_YI,f_ZI);
g2 = surf(g_XI,g_YI,g_ZI);

set(g1,'FaceColor',[1 0 0],'FaceAlpha',0.9);
set(g2,'FaceColor',[0 1 0],'FaceAlpha',0.9);

%alpha(g3,0.05)
hold off


end


function [fxy_coordinates] = getCoordinates_CP(fxy_matrix_Brn)
%% Get the (x,y,z) coordinates of each of the control points
%
% % Inputs
%
% fxy_matrix_Brn : (Matrix) Coefficients of f(x,y)

% get dimensions of f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy_matrix_Brn);

% Initialise a matrix to store all of the coordinates
num_coordinates = (m1+1) * (m2+1);
fxy_coordinates = zeros(num_coordinates, 3);


count = 1;
for i = 0:1:m1
    for j = 0:1:m2
        % get x coordinate of control point
        x_pos = i/(m1);
        % Get y coordinate of control point
        y_pos = j/(m2);
        % Get z coordinate of control point
        z_pos = fxy_matrix_Brn(i+1,j+1);
        
        % Create new coordinate
        new_coord = [x_pos y_pos z_pos];
        
        % Add the coordinates to the dataset
        fxy_coordinates(count,:) = new_coord;
        
        count = count + 1;
    end
end


end

function fxy_coordinates_Curve = get_Coordinates_Curve(fxy_matrix_Brn)
% Given a set of coefficients for a surface f(x,y) obtain a set of (x,y,z)
% values for points on the curve.
%
% % Inputs
%
% fxy_matrix_Brn : (Matrix) Coefficients of f(x,y)


[m1, m2] = GetDegree_Bivariate(fxy_matrix_Brn);

% Initialise empty vector
fxy_coordinates_Curve = [];

% % Evaluate the Bernstein surface over the interval
for i = 0:0.1:m1
    x_inc = i/m1;
    for j = 0:0.1:m2
        y_inc = j/m2;
        
        % Obtain the z component
        z_Brn_exact = Evaluate_BernsteinPoly(x_inc,y_inc,...
            fxy_matrix_Brn);
        
        % Add the (x,y,z) to the set of coordinates
        fxy_coordinates_Curve = [fxy_coordinates_Curve; x_inc y_inc z_Brn_exact];
        
    end
end

end
