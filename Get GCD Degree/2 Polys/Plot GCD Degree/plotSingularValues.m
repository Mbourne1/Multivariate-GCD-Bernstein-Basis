function [] = plotSingularValues(arrSingularValues, myLimits_t1, myLimits_t2, limits_t1, limits_t2)
%
% Plot all of the singular values of S_{k1,k2} for k1 =
% lowerlim,...,upperlim and k2 = lowerlim,...,upperlim
%
% % Inputs
%
% arrSingularValues : 2 Dimensional array stores vectors containing all
% singular values of S_{k1,k2}.
%
% myLimits_t1 :
%
% myLimits_t2 :
%
% limits_t1 :
%
% limits_t2 :

% Get upper and lower limits for t1 and t2
lowerLimit_t1 = myLimits_t1(1);
upperLimit_t1 = myLimits_t1(2);

lowerLimit_t2 = myLimits_t2(1);
upperLimit_t2 = myLimits_t2(2);



% Get number of subresultants whose singular values are to be plotted
nSubresultants_t1 = upperLimit_t1 - lowerLimit_t1 +1;
nSubresultants_t2 = upperLimit_t2 - lowerLimit_t2 +1;

global SETTINGS
figure_name = sprintf('Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on

for i1 = 1:1:nSubresultants_t1    
    for i2 = 1:1:nSubresultants_t2

        % Get k1, k2 subscript of the Sylvester matrix in array location
        % i1,i2
        k1 = lowerLimit_t1 + (i1-1);
        k2 = lowerLimit_t2 + (i2-1);
        
        % Get vector of singular values
        vSingularValues = arrSingularValues{i1, i2};
        
        v_k1 = k1.*ones(length(vSingularValues));
        v_k2 = k2.*ones(length(vSingularValues));
        
        plot3(v_k1, v_k2, log10(vSingularValues),'*');
        
    end
end
grid on
hold off

end
