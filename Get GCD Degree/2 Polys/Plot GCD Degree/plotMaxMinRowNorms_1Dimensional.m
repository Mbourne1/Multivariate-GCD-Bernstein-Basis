function [] = plotMaxMinRowNorms_degreeRelative_1Dimensional(vMaxRowNorm, vMinRowNorm, my_limits_t1)
%
% % Inputs
%
% vMaxRowNorm : 
%
% vMinRowNorm : 
%
% my_limits_t1 :
%

global SETTINGS

% Get upper and lower bounds
lowerLimit_t1 = my_limits_t1(1);
upperLimit_t1 = my_limits_t1(2);
vec_x = lowerLimit_t1 : 1 : upperLimit_t1;
temp_vec = (vMinRowNorm) ./ (vMaxRowNorm);
figure_name = sprintf('Plot Max:Min Row Norms of %s',SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name', figure_name)
hold on
mesh(vec_x, log10(temp_vec') );
xlabel('k_{1}')
ylabel('Min/Max Row Norm of R1')
grid on
hold off


end

