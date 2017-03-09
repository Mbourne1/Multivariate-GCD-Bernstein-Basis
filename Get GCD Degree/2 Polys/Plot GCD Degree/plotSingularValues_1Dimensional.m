function [] = plotSingularValues_1Dimensional(arrSingularValues, myLimits_t1, limits_t1)
%
%
% % Input
%
% arrSingularValues : (1D Array of Vectors) : Vectors containing singular
% values of the Sylvester subresultant matrices
%
% my_limits_t1 : [(Int) (Int)] Pair of integers to denote the lower and
% upper bound on the possible values of the degree of the GCD.
%

global SETTINGS

% Get lower and upper bound
myLowerLimit_t1 = myLimits_t1(1);
myUpperLimit_t1 = myLimits_t1(2);

lowerLimit_t1 = limits_t1(1);
upperLimit_t1 = limits_t1(2);

% Get number of Sylvester subresultants
nSubresultants_t1 = myUpperLimit_t1 - myLowerLimit_t1 +1;

% Plot Figure
figure_name = sprintf('Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on

for i1 = 1:1:nSubresultants_t1
   

        k1 = myLowerLimit_t1 + (i1-1);

        % Get vector of singular values 
        vSingularValues = arrSingularValues{i1};
        v_k1 = k1.*ones(length(vSingularValues));
        plot(v_k1, log10(vSingularValues),'*');
        

end
vline(lowerLimit_t1);
vline(upperLimit_t1);
grid on
hold off

end
