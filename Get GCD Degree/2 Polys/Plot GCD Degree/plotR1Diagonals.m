function [] = plotR1Diagonals(arr_R1, limits_k1, limits_k2, limits_t1, ...
    limits_t2)
% Plot the diagonals of the matrices R_{1,k1,k2} obtained by QR decomposition
% of the subresultant matrices S_{k1,k2}
%
% % Inputs
%
% arr_R1 : (Array of Matrices) Array of matrices R_{1,k1,k2} obtained by QR
% decomposition of S_{k1,k2}.
%
% limits_k1 : (Int Int) Range of values of k_{1}
%
% limits_k2 : (Int Int) Range of values of k_{2}
%
% limits_t1 : (Int Int) Range in which t_{1} (The degree of the GCD with
% respect to x) is known to lie.
%
% limits_t2 : (Int Int) Range in which t_{2} (The degree of the GCD with
% respect to y) is known to lie.

% Get upper and lower limits
lowerLimit_k1 = limits_k1(1);
upperLimit_k1 = limits_k1(2);

lowerLimit_k2 = limits_k2(1);
upperLimit_k2 = limits_k2(2);

% Get number of Subresultant matrices
nSubresultants_t1 = upperLimit_k1 - lowerLimit_k1 +1 ;
nSubresultants_t2 = upperLimit_k2 - lowerLimit_k2 +1 ;

global SETTINGS

% Figure Head
figure_name = sprintf('Diagonals of R1 from QR decomp of %s',SETTINGS.SYLVESTER_MATRIX_VARIANT);
figure('name', figure_name)
hold on

% For each subresultant matrix
for i1 = 1:1:nSubresultants_t1
    for i2 = 1:1:nSubresultants_t2
        
        k1 = lowerLimit_k1 + (i1 -1);
        k2 = lowerLimit_k2 + (i2 -1);
        
        % Get diagonal of R_{1, k_{1}, k_{2}}
        temp_vec = diag(arr_R1{i1, i2});
        
        v_k1 = k1 .* ones(length(temp_vec), 1);
        v_k2 = k2 .* ones(length(temp_vec), 1);
        
        plot3(v_k1, v_k2, log10(temp_vec), '*');
        
    end
end

% Figure Labels
xlabel('k_{1}')
ylabel('k_{2}')
zlabel('log_{10} Row Diagonals of R_{1}')

% Figure Layout
grid on
hold off

end