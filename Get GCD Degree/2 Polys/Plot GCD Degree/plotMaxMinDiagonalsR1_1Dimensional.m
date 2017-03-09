function [] = plotMaxMinDiagonalsR1_1Dimensional(vec_MaxDiagonals_R1, vec_MinDiagonals_R1, myLimits_t1, limits_t1)
%
% % Inputs
%
% mat_MaxDiagonals_R1 : (Matrix)
%
% mat_MinDiagonals_R1 : (Matrix)
%
% myLimits_t1 : ([Int Int])
%
% limits_t1 :


% Get upper and lower limits of k1 and k2
myLowerLimit_k1 = myLimits_t1(1);
myUpperLimit_k1 = myLimits_t1(2);

%
lowerLimit_k1 = limits_t1(1);
upperLimit_k1 = limits_t1(1);

%
v_k1 = myLowerLimit_k1:1:myUpperLimit_k1;


temp_vec = log10(vec_MinDiagonals_R1./vec_MaxDiagonals_R1);

figure_name = 'Plot Max:Min Values';
figure('name',figure_name)
hold on
plot(v_k1, temp_vec,'-s');
xlim([myLowerLimit_k1, myUpperLimit_k1]);
xlabel('k_{1}');
ylabel('log_{10} Min/Max Diagonals of R_{1}')
title(figure_name);
vline(lowerLimit_k1);
vline(upperLimit_k1);
grid on
hold off

end