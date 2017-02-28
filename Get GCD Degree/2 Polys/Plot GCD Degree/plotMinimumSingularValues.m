function [] = plotMinimumSingularValues_degreeRelative(mat_MinimumSingularValues, my_limits_t1, my_limits_t2)
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
lowerLimit_t2 = my_limits_t2(1);
upperLimit_t2 = my_limits_t2(2);


i1 = lowerLimit_t1:1:upperLimit_t1;
i2 = lowerLimit_t2:1:upperLimit_t2;
[X,Y] = meshgrid(i1,i2);

figure_name = sprintf('Minimum Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on
title('Minimum Singular Values');
xlabel('k_{1}');
ylabel('k_{2}');
grid on
surf(X,Y,log10(mat_MinimumSingularValues)');
hold off


end


