function [] = plotMaxMinDiagonalsR1_1Dimensional(vMaxDiagonals_R1, ...
    vMinDiagonals_R1, limits_k1, limits_t1)
%
% % Inputs
%
% vMaxDiagonals_R1 : (Matrix)
%
% vMinDiagonals_R1 : (Matrix)
%
% limits_k1 : (Int Int) Upper and lower bound for k_{1}
%
% limits_t1 : (Int Int) Upper and lower bound for t_{1}


% Get upper and lower limits of k_{1}
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

% Get upper and lower limit of t_{1}
lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(1);

% Get vector of all possible k_{1}
v_k1 = lowerLimit_k1:1:upperLimit_k1;

% Get rank revealing metric
metric = log10(vMinDiagonals_R1./vMaxDiagonals_R1);

% Plot
figure_name = 'Plot Max:Min Values';
figure('name',figure_name)
hold on
plot(v_k1, metric,'-s');
xlim([lowerLimit_k1, upperLimit_k1]);
xlabel('k_{1}');
ylabel('log_{10} Min/Max Diagonals of R_{1}')
title(figure_name);

% Plot upper and lower bound for t_{1}
vline(lowerLimit_t1);
vline(upperLimit_t1);

grid on
hold off

end