function [] = plotSingularValues(arr_SingularValues, my_limits_t1, my_limits_t2)

lowerLimit_t1 = my_limits_t1(1);
upperLimit_t1 = my_limits_t1(2);
lowerLimit_t2 = my_limits_t2(1);
upperLimit_t2 = my_limits_t2(2);

nSubresultants_t1 = upperLimit_t1 - lowerLimit_t1 +1;
nSubresultants_t2 = upperLimit_t2 - lowerLimit_t2 +1;

global SETTINGS

figure_name = sprintf('Singular Values of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
hold on

for i1 = 1:1:nSubresultants_t1
    
    for i2 = 1:1:nSubresultants_t2

        k1 = lowerLimit_t1 + (i1-1);
        k2 = lowerLimit_t2 + (i2-1);
        
        % Get vector of singular values
        vSingularValues = arr_SingularValues{i1, i2};
        
        v_k1 = k1.*ones(length(vSingularValues));
        v_k2 = k2.*ones(length(vSingularValues));
        
        plot3(v_k1, v_k2, log10(vSingularValues),'*');
        
    end
end
grid on
hold off

end
