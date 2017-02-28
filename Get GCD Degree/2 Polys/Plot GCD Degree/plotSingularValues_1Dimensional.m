function [] = plotSingularValues_1Dimensional(arrSingularValues, my_limits_t1)
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
lowerLimit_t1 = my_limits_t1(1);
upperLimit_t1 = my_limits_t1(2);

% Get number of Sylvester subresultants
nSubresultants_t1 = upperLimit_t1 - lowerLimit_t1 +1;

% Plot Figure
figure_name = sprintf('Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on

for i1 = 1:1:nSubresultants_t1
   

        k1 = lowerLimit_t1 + (i1-1);

        % Get vector of singular values 
        vSingularValues = arrSingularValues{i1};
        v_k1 = k1.*ones(length(vSingularValues));
        plot(v_k1, log10(vSingularValues),'*');
        

end
grid on
hold off

end
