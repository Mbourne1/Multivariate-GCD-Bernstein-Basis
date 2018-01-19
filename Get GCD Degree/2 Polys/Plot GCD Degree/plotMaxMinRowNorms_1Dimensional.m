function [] = plotMaxMinRowNorms_1Dimensional(vMaxRowNorm, vMinRowNorm,...
    limits_k1, limits_t1)
%
% % Inputs
%
% vMaxRowNorm : (Vector) 
%
% vMinRowNorm : (Vector)
%
% limits_k1 : (Int Int)
%
% limits_t1 : (Int Int)
%


global SETTINGS

% Get upper and lower bounds
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

% Get upper and lower limits
lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);

% Get vector x
vec_x = lowerLimit_k1 : 1 : upperLimit_k1;

temp_vec = (vMinRowNorm) ./ (vMaxRowNorm);
figure_name = sprintf('Plot Max:Min Row Norms of %s',SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name', figure_name)
hold on

mesh(vec_x, log10(temp_vec') );
xlabel('k_{1}')
ylabel('Min/Max Row Norm of R1')
vline(lowerLimit_t1);
vline(upperLimit_t1);
grid on

hold off


end

