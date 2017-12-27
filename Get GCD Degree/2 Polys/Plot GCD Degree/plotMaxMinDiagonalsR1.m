function [] = plotMaxMinDiagonalsR1(mat_MaxDiagonals_R1, ...
    mat_MinDiagonals_R1, limits_k1, limits_k2, limits_t1, limits_t2)
%
% % Inputs
%
% mat_MaxDiagonals_R1 : (Matrix)
%
% mat_MinDiagonals_R1 : (Matrix)
%
% myLimits_t1 : ([Int Int])
%
% myLimits_t2 : ([Int Int])
%
% limits_t1 : (Int Int)
% 
% limits_t2 : (Int Int)

% Get upper and lower limits of k1 and k2
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

v_k1 = lowerLimit_k1 : 1 : upperLimit_k1;
v_k2 = lowerLimit_k2 : 1 : upperLimit_k2;

[X,Y] = meshgrid(v_k1,v_k2);

mat = log10(mat_MinDiagonals_R1./mat_MaxDiagonals_R1);

figure_name = 'Plot Max:Min Values';
figure('name',figure_name)
hold on

mesh(X,Y,mat');
xlim([lowerLimit_k1, upperLimit_k1]);
ylim([lowerLimit_k2, upperLimit_k2]);
xlabel('k_{1}');
ylabel('k_{2}');
zlabel('log_{10} Min/Max Diagonals of R_{1}')
title(figure_name);
grid on
hold off

end