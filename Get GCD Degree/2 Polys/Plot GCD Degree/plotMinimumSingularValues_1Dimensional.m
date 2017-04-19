function [] = plotMinimumSingularValues_1Dimensional(vMinimumSingularValues, limits_k1, limits_t1)
%
% % Inputs
%
% vMinimumSingularValues : (Vector)
%
% limits_k1 : (Int Int)
%
% limits_t1 : (Int Int)

global SETTINGS

% Get upper and lower limits for k_{1}
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

% Get upper and lower limits for t_{1}
lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);

%
vec_x = lowerLimit_k1 : 1 : upperLimit_k1;

figure_name = sprintf('Minimum Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on
title('Minimum Singular Values');

vline(lowerLimit_t1);
vline(upperLimit_t1);

grid on
plot(vec_x,log10(vMinimumSingularValues)');
hold off


end


