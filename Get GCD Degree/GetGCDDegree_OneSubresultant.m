function t = GetGCDDegree_OneSubresultant(vMetric, myLimits)
% Given the vector of values from either minimum singular values or max:min
% R diagonals.
% Get the rank, where only one subresultant exists.
%
% % Inputs
%
% vMetric : (Vector)
%
% myLimits : [(Int) (Int)]
%
% % Outputs
%
% t : Degree of Greatest common divisor

myLowerLimit = myLimits(0);
myUpperLimit = myLimits(1);

global SETTINGS

% Get the calling function
[St,~] = dbstack();
calling_function = St(2).name;

% Only one subresultant
fprintf([calling_function ' : ' mfilename ' : ' 'Only one subresultant exists. \n'])


if(SETTINGS.PLOT_GRAPHS)
        figure_name = sprintf([calling_function ' : Singular values of S_{1} \n']);
        figure('name',figure_name)
        hold on
        title('Singular values of S_{1}')
        plot(log10(vMetric))
        hold off
    
end

[deltaSingularValues,~] = Analysis(vMetric);

% If the change is smaller than the predefined threshold value, then plot
% is considered 'flat'.
if deltaSingularValues < SETTINGS.THRESHOLD
    
    % The subresultant is of full rank, in which case t = 0
    t = 0;
    fprintf([calling_function ' : ' 'The only Subresultant S_{1} appears to be of NonSingular. \n']);
    return
    
else % val > threshold
    
    % The subresultant S_{1} is rank deficient, in which case t = 1
    t = myUpperLimit;
    fprintf([calling_function ' : ' 'The only Subresultant S_{1} appears to be Singular \n']);
    return
    
end
end
