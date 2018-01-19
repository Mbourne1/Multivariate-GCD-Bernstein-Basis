function [] = plotMaxMinRowNorms(mat_MaxRowNorm, mat_MinRowNorm, ...
    limits_k1, limits_k2, limits_t1, limits_t2)
%
% % Inputs
%
% mat_MaxRowNorm : (Matrix)
%
% mat_MinRowNorm : (Matrix)
%
% limits_k1 : (Int Int) 
%
% limits_k2 : (Int Int) 
%
% limits_t1 : (Int Int) 
%
% limits_t2 : (Int Int)

global SETTINGS

% Get upper and lower bounds
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

lowerLimit_k2 = limits_k1(1);
upperLimit_k2 = limits_k2(2);


vec_x = lowerLimit_k1 : 1 : upperLimit_k1;
vec_y = lowerLimit_k2 : 1 : upperLimit_k2;

[X, Y] = meshgrid(vec_x, vec_y);

mat = (mat_MinRowNorm) ./ (mat_MaxRowNorm);

figure_name = sprintf('Plot Max:Min Row Norms of %s',SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name', figure_name)
hold on
mesh(X, Y, log10(mat') );
xlabel('k_{1}')
ylabel('k_{2}')
zlabel('Min/Max Row Norm of R1')
grid on
hold off


end

