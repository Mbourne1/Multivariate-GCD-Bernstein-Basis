function [] = plotMinimumSingularValues_1Dimensional( ...
    vMinimumSingularValues, limits_k, limits_t, ...
    txtTitle, txtXLabel, txtYLabel)
%
% % Inputs
%
% vMinimumSingularValues : (Vector) Vector of minimum singular values of
% the set of subresultant matrices S_{k}(f,g)
%
% limits_k1 : (Int Int) Range of k values
%
% limits_t1 : (Int Int) Set of k in which the degree of the GCD (t) lies.
%
% txtTitle : (String)
%
% txtXLabel : (String)
%
% txtYLabel : (String)


global SETTINGS

% Get upper and lower limits for k_{1}
lowerLimit_k1 = limits_k(1);
upperLimit_k1 = limits_k(2);

% Get upper and lower limits for t_{1}
lowerLimit_t1 = limits_t(1);
upperLimit_t1 = limits_t(2);

% Initialise vector x
vec_x = lowerLimit_k1 : 1 : upperLimit_k1;

% Figure Head
figure_name = sprintf('Minimum Singular Values of %s', SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name',figure_name)
hold on
title('Minimum Singular Values');

%vline(lowerLimit_t1);
%vline(upperLimit_t1);

title(txtTitle);
xlabel(txtXLabel, 'Interpreter', 'latex');
ylabel(txtYLabel, 'Interpreter', 'latex');

plot(vec_x, log10(vMinimumSingularValues)', '-s');
hold off

% Figure Layout
grid on

end


