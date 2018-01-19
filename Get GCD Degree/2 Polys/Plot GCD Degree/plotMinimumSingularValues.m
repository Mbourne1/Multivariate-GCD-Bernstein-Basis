function [] = plotMinimumSingularValues(matMinimumSingularValues, ...
    limits_k1, limits_k2, limits_t1, limits_t2, rank_range)
%
% % Inputs
%
% matMinimumSingularValues : (Matrix) Stores the minimum singular values of
% each Sylvester subresultant matrix S_{k1,k2} for k1 =
% lowerlim,...,upperlim and k2 = lowerlim,...,upperlim
%
% limits_k1 : [Int Int]
%
% limits_k2 : [Int Int]
%
% limits_t1 : [Int Int]
%
% limits_t2 : [Int Int]
%
% rank_range : [Float Float]

global SETTINGS

% Get upper and lower limit
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

%
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

%
v_i1 = lowerLimit_k1 : 1 : upperLimit_k1;
v_i2 = lowerLimit_k2 : 1 : upperLimit_k2;

[X,Y] = meshgrid(v_i1,v_i2);

figure_name = sprintf('Minimum Singular Values of %s', SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name',figure_name)
hold on
%title('Minimum Singular Values');

try
    surf(X, Y, log10(matMinimumSingularValues)');
    
    % Labels
    xlabel('$k_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$k_{2}$', 'Interpreter', 'latex', 'FontSize', 20);
    zlabel('$\log_{10} \left( \dot{\sigma}_{k_{1}, k_{2}} \right)$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Set view angle
    az = -30;
    el = 20;
    view(az, el);
    
    
    % Display
    grid on
    box on
    
    % Set location of window and size
    m_left = 100;
    m_bottom = 100;
    m_width = 600;
    m_height = 600;
    
    set(gcf, 'Position', [m_left, m_bottom, m_width, m_height]);
    
    % Position of figure within window
    myplot = gca;
    myval_side = 0.12;
    myval_base = 0.10;
    set(myplot, 'Position', [ myval_side myval_base 0.98 - myval_side 0.98 - myval_base])
    
    myText = sprintf("emin : %e",  SETTINGS.EMIN);
    
    dim = [.1 .65 .3 .3];
    %str = 'Straight Line Plot from 1 to 10';
    annotation('textbox',dim,'String',myText,'FitBoxToText','on');
    
   
    
    hold off
catch
    
    fprintf('Surface Not Plotted \n')
end

end


