function [] = plotMaxMinRowNorms(mat_MaxRowNorm, mat_MinRowNorm, myLimits_t1, myLimits_t2, limits_t1, limits_t2)
%
% % Inputs
%
% mat_MaxRowNorm : 
%
% mat_MinRowNorm : 
%
% myLimits_t1 :
%
% myLimits_t2 :
%
% limits_t1 :
%
% limits_t2 :

global SETTINGS

% Get upper and lower bounds
lowerLimit_t1 = myLimits_t1(1);
upperLimit_t1 = myLimits_t1(2);
lowerLimit_t2 = myLimits_t2(1);
upperLimit_t2 = myLimits_t2(2);

vec_x = lowerLimit_t1 : 1 : upperLimit_t1;
vec_y = lowerLimit_t2 : 1 : upperLimit_t2;

[X, Y] = meshgrid(vec_x, vec_y);

mat = (mat_MinRowNorm) ./ (mat_MaxRowNorm);

figure_name = sprintf('Plot Max:Min Row Norms of %s',SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name', figure_name)
hold on
mesh(X, Y, log10(mat') );
xlabel('k_{1}')
ylabel('k_{2}')
zlabel('Min/Max Row Norm of R1')
grid on
hold off


end

