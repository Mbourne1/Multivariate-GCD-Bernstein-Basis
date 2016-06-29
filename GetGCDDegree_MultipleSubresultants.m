function [t] = GetGCDDegree_MultipleSubresultants(vMinimumSingularValues,deg_limits)

global SETTINGS

% Get the name of the function which called this function
[St,~] = dbstack();
calling_function = St(2).name;

% Get maximum change in the minimum singular values
[maxDelta,index] = Analysis(vMinimumSingularValues);

switch SETTINGS.PLOT_GRAPHS
    case 'y'
        plot(log10(vMinimumSingularValues));
    case 'n'
end

% Get upper and lower bounds
lower_lim = deg_limits(1);
upper_lim = deg_limits(2);

% check if the maximum change is significant
fprintf([mfilename ' : ' sprintf('Threshold :  %2.4f \n', SETTINGS.THRESHOLD)]);
fprintf([mfilename ' : ' sprintf('Max change : %2.4f \n', maxDelta)]);


if abs(maxDelta) < SETTINGS.THRESHOLD
    % Change in Singular values is not significant so check if all
    % subresultants are rank deficient or full rank
    
    % Get the average of the minimum singular values
    avg = mean(vMinimumSingularValues);
    
    if avg < SETTINGS.THRESHOLD_RANK
       % All Minimum singular values are below threshold so, all 
       % subresultants are rank deficient. deg(GCD) = 0
       fprintf([calling_function ' : ' mfilename ' : ' 'Polynomails are coprime\n' ])
       t = 0;
    else 
        % All minimum singular values are above threshold so all
        % subresultants are full rank. deg(GCD) = min(m,n)
       fprintf([calling_function ' : ' mfilename ' : ' 'All Subresultants are rank deficient \n' ])
       t = upper_lim;
    end

else
    % change is significant
    fprintf([mfilename ' : ' 'Significant Change' ]);
    t = lower_lim + index - 1;
    fprintf(': %i \n',t);
    
end

end