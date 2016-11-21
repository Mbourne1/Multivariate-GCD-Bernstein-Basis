function [] = Plot_Parametric(parametric_arr)
% Inputs 

% parametric_arr : contains 'n' sets of data points for 'n' parametrically
% defined curves.

% Get number of parametric curves in the parametric curve array
[r,n] = size(parametric_arr);



figure('name','Parametrically Defined Curves')
hold on
% For every parametrically defined curve
for i = 1:1:n
    
    % Get the ith parametric curve
    parametric = parametric_arr{i};
    
    % Get the equation in x(t)
    X = parametric(1,:);
    % Get the equation  y(t)
    Y = parametric(2,:);
    
    name = sprintf('Curve C_{%i}',i)
    plot(X,Y,'DisplayName',name)
    
end
grid on
legend('show')
hold off




end