function [] = PlotR1RowNorms_degreeRelative_1Dimensional(arr_R1_RowNorms, my_limits_t1)
%             
% % Inputs
%
% arr_R1_RowNorms : (Array)
%
% my_limits_t1 : ([Int Int])



global SETTINGS

% Get limits
lowerLimit_t1 = my_limits_t1(1);
upperLimit_t1 = my_limits_t1(2);

% Get number of subresultants
nSubresultants_t1 = upperLimit_t1 - lowerLimit_t1 + 1;

% Plot
figure_name = sprintf('R1 Row Norms of QR decomposition of %s', SETTINGS.SYLVESTER_BUILD_METHOD);
figure('name',figure_name)
title(figure_name)
hold on

for i1 = 1:1:nSubresultants_t1
        
        k1 = lowerLimit_t1 + (i1 - 1);
  
        temp_vec = arr_R1_RowNorms{i1,i2};
        
        vec_k1 = k1.* ones(length(temp_vec),1);
       
        
        plot2(vec_k1, log10(temp_vec),'*');
        

end
xlabel('k')
zlabel('Row Norms of R_{1}')
grid on
hold off

end