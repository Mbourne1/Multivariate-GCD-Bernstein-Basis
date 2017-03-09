function [] = plotMinimumSingularValues_1Dimensional(vMinimumSingularValues, myLimits_t1, limits_t1)
%
% % Inputs
%
% mat_MinimumSingularValues : 
%
% myLimits_t1
%
% limits_t1

global SETTINGS

%
myLowerLimit_t1 = myLimits_t1(1);
myUpperLimit_t1 = myLimits_t1(2);

%
lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);

%
vec_x = myLowerLimit_t1 : 1 : myUpperLimit_t1;

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


