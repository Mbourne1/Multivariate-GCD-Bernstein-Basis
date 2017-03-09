function [] = plotMaxMinRowNorms_1Dimensional(vMaxRowNorm, vMinRowNorm, myLimits_t1, limits_t1)
%
% % Inputs
%
% vMaxRowNorm : 
%
% vMinRowNorm : 
%
% myLimits_t1 :
%
% limits_t1 :

global SETTINGS

% Get upper and lower bounds
myLowerLimit_t1 = myLimits_t1(1);
myUpperLimit_t1 = myLimits_t1(2);

% Get upper and lower limits
lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);

% Get vector x
vec_x = myLowerLimit_t1 : 1 : myUpperLimit_t1;

temp_vec = (vMinRowNorm) ./ (vMaxRowNorm);
figure_name = sprintf('Plot Max:Min Row Norms of %s',SETTINGS.SYLVESTER_BUILD_METHOD);
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

