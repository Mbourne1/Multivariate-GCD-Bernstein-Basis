function [] = PlotImplicitBezierSurface(fxy_matrix_Brn)
% Plot the Bezier surfaces whose control points are uniformly distributed
% over the x and y axes.
% Given a bivariate polynomial, plot the surfaces

% % Inputs

% fxy_matrix_Brn:

%               B_{0}(y)  B_{1}(y) ....
%                _______ _______
%   B_{0}(x)    |_______|_______| .....
%   B_{1}(x)    |_______|_______|
%    ...        |_______|_______|
%   B_{m1}(x)   |_______|_______|

%% Plot 1

% % Plot the control points of the surface in bernstein basis
[m1,m2] = GetDegree_Bivariate(fxy_matrix_Brn);

% Get number of control points
num_cp = (m1+1) * (m2+1);

% Store x y and z of each control point
data = zeros(num_cp,3);

% Set counter to one
count = 1;
for i = 0:1:m1
    for j = 0:1:m2
        % get x coordinate of control point
        x_pos = i/(m1);
        % Get y coordinate of control point
        y_pos = j/(m2);
        % Get z coordinate of control point
        z_pos = fxy_matrix_Brn(i+1,j+1);
        % Add the coordinates to the dataset
        data(count,:) = [x_pos y_pos z_pos];
        z(i+1,j+1) = z_pos;
        count = count + 1;
    end
end

% Plot the control points in a three dimensional scatter
figure('name','Bernstein Control Points', 'NumberTitle','off')
hold on
scatter3(data(:,1), data(:,2), data(:,3),'Filled')
hold off

%Get the increment size in x
inc_x = 1/m1;
inc_y = 1/m2;

% Generate the mesh grid
x = 0:inc_x:1;
y = 0:inc_y:1;
[x,y] = meshgrid(x,y);

figure('name','Surface Plot')
hold on
s1 = surf(x,y,z');
xlabel('t_{1}')
ylabel('t_{2}')
alpha(s1,0.5)
xlabel('t_{1}')
ylabel('t_{2}')
hold off
%% Plot 3
% % Plot the Bézier surface by evaluating at a series of points

data_Brn = [];

% % Evaluate the Bernstein surface over the interval
for i = 0:0.1:m1
    x_inc = i/m1;
    for j = 0:0.1:m2
        y_inc = j/m2;
        
        z_Brn_exact = Evaluate_BernsteinPoly(x_inc,y_inc,...
            fxy_matrix_Brn);
        
        z_surf_data = z_Brn_exact;
        
        data_Brn = [data_Brn; x_inc y_inc z_Brn_exact];
        
    end
end


% % plot the exact Bernstein data and the noisy bernstein data
figure('name','Bernstein Surfaces','NumberTitle','off')
title('Bernstein Surface')
hold on

% Plot exact data
scatter3(data_Brn(:,1),data_Brn(:,2),data_Brn(:,3),'Filled')

hold off