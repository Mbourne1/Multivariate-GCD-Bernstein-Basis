function [t] = GetGCDDegree_MultipleSubresultants(vMetric, myLimits_t)
%
% % Inputs
%
% vMetric : (Vector) Contains value which identifies whether corresponding
% Sylvester matrix is full rank or rank deficient. Low values indicate rank
% deficiency, high values indicate full rank.
%
% myLimits_t : [(Int) (Int)] 
%
% % Outputs
%
% t : Degree of greatest common divisor


global SETTINGS


% Get the name of the function which called this function
[St,~] = dbstack();
calling_function = St(2).name;

% Get maximum change in the vector of rank revaling metric values
[maxDelta, idxMaxChange] = Analysis(vMetric);

% Get upper and lower bounds
myLowerLimit = myLimits_t(1);
myUpperLimit = myLimits_t(2);

% check if the maximum change is significant
fprintf([mfilename ' : ' sprintf('Threshold :  %2.4f \n', SETTINGS.THRESHOLD)]);
fprintf([mfilename ' : ' sprintf('Max change : %2.4f \n', maxDelta)]);


if abs(maxDelta) < SETTINGS.THRESHOLD
    % Change in Singular values is not significant so check if all
    % subresultants are rank deficient or full rank
    
    % Get the average of the minimum singular values
    avg = mean(vMetric);
    
    if avg < SETTINGS.THRESHOLD_RANK
        
        % All Minimum singular values are below threshold so, all
        % subresultants are rank deficient. deg(GCD) = 0
        fprintf([calling_function ' : ' mfilename ' : ' 'Polynomails are coprime\n' ])
        t = 0;
        
    else
        
        % All minimum singular values are above threshold so all
        % subresultants are full rank. deg(GCD) = min(m,n)
        fprintf([calling_function ' : ' mfilename ' : ' 'All Subresultants are rank deficient \n' ])
        t = myUpperLimit;
        
    end
    
else
    
    % Delta value is significant
    fprintf([mfilename ' : ' 'Significant Change detected \n' ]);
    t = myLowerLimit + (idxMaxChange - 1);
   
    
end

end