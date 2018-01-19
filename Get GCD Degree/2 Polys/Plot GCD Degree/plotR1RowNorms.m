function [] = plotR1RowNorms(arr_R1_RowNorms, limits_k1, limits_k2, limits_t1, limits_t2)
%
% % Inputs
%
% arr_R1_RowNorms : (Array)
%
% limits_k1 : (Int Int) 
%
% limits_k2 : (Int Int)
%
% limits_t1 : (Int Int)
%
% limits_t2 : (Int Int)

% Note that limits_t1 and limits_t2 are not currently used, but I want to
% include lines on the plot to show where these values are.

global SETTINGS

% Get limits
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);
lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

% Get number of subresultants
nSubresultants_t1 = upperLimit_k1 - lowerLimit_k1 + 1;
nSubresultants_t2 = upperLimit_k2 - lowerLimit_k2 + 1;

% Plot
figure_name = sprintf('R1 Row Norms of QR decomposition of %s', SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name',figure_name)
title(figure_name)
hold on

for i1 = 1:1:nSubresultants_t1
    for i2 = 1:1:nSubresultants_t2
        
        k1 = lowerLimit_k1 + (i1 - 1);
        k2 = lowerLimit_k2 + (i2 - 1);
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