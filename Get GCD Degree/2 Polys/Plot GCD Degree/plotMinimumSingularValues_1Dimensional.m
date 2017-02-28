function [] = plotMinimumSingularValues_degreeRelative_1Dimensional(vMinimumSingularValues, my_limits_t1)
%
% % Inputs
%
% mat_MinimumSingularValues : 
%
% my_limits_t1
%
% my_limits_t2

global SETTINGS

lowerLimit_t1 = my_limits_t1(1);
upperLimit_t1 = my_limits_t1(2);

vec_x = lowerLimit_t1:1:upperLimit_t1;



figure_name = sprintf('Minimum Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on
title('Minimum Singular Values');

grid on
plot(vec_x,log10(vMinimumSingularValues)');
hold off


end


