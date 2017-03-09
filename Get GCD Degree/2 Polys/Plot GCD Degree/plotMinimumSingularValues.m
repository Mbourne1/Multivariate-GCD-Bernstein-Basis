function [] = plotMinimumSingularValues(matMinimumSingularValues, myLimits_t1, myLimits_t2, limits_t1, limits_t2)
%
% % Inputs
%
% matMinimumSingularValues : (Matrix) Stores the minimum singular values of
% each Sylvester subresultant matrix S_{k1,k2} for k1 =
% lowerlim,...,upperlim and k2 = lowerlim,...,upperlim
%
% my_limits_t1 :
%
% my_limits_t2 :
%
% limits_t1 :
%
% limits_t2

global SETTINGS

% Get upper and lower limit
lowerLimit_t1 = myLimits_t1(1);
upperLimit_t1 = myLimits_t1(2);
lowerLimit_t2 = myLimits_t2(1);
upperLimit_t2 = myLimits_t2(2);


v_i1 = lowerLimit_t1 : 1 : upperLimit_t1;
v_i2 = lowerLimit_t2 : 1 : upperLimit_t2;

[X,Y] = meshgrid(v_i1,v_i2);

figure_name = sprintf('Minimum Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on
title('Minimum Singular Values');
xlabel('k_{1}');
ylabel('k_{2}');
grid on
surf(X,Y,log10(matMinimumSingularValues)');
hold off


end


