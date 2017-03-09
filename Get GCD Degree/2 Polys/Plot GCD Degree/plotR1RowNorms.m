function [] = plotR1RowNorms(arr_R1_RowNorms, myLimits_t1, myLimits_t2, limits_t1, limits_t2)
%
% % Inputs
%
% arr_R1_RowNorms : (Array)
%
% my_limits_t1 : ([Int Int])
%
% my_limits_t2 : ([Int Int])


global SETTINGS

% Get limits
lowerLimit_t1 = myLimits_t1(1);
upperLimit_t1 = myLimits_t1(2);
lowerLimit_t2 = myLimits_t2(1);
upperLimit_t2 = myLimits_t2(2);

% Get number of subresultants
nSubresultants_t1 = upperLimit_t1 - lowerLimit_t1 + 1;
nSubresultants_t2 = upperLimit_t2 - lowerLimit_t2 + 1;

% Plot
figure_name = sprintf('R1 Row Norms of QR decomposition of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
title(figure_name)
hold on

for i1 = 1:1:nSubresultants_t1
    for i2 = 1:1:nSubresultants_t2
        
        k1 = lowerLimit_t1 + (i1 - 1);
        k2 = lowerLimit_t2 + (i2 - 1);
        temp_vec = arr_R1_RowNorms{i1,i2};
        
        vec_k1 = k1.* ones(length(temp_vec),1);
        vec_k2 = k2.* ones(length(temp_vec),1);
        
        plot3(vec_k1, vec_k2, log10(temp_vec),'*');
        
    end
end
xlabel('k_{1}')
ylabel('k_{2}')
zlabel('Row Norms of R_{1}')
grid on
hold off

end