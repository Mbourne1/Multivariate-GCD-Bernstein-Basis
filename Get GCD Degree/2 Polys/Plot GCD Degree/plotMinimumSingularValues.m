function [] = plotMinimumSingularValues(matMinimumSingularValues, limits_k1, limits_k2, limits_t1, limits_t2, rank_range)
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

figure_name = sprintf('Minimum Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on
title('Minimum Singular Values');
xlabel('k_{1}');
ylabel('k_{2}');
grid on
surf(X, Y, log10(matMinimumSingularValues)');
hold off


end


