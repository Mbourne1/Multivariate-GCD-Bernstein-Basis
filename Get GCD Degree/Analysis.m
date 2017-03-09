function [max_Delta, idxMaxChange] = Analysis(vMetric)
% 
% % Inputs
%
% vMetric : (Vector) Values from Rank revealing metric. Low values suggest
% low rank, high values suggest full rank.
%
% % Outputs
%
% max_Delta : Maximum change
%
% idxMaxChange : Location of maximum change in the vector


    vMetric = sanitize(vMetric);
    
    % Get the change in the ratios from one subresultant to the next.
    vDelta = abs(diff(log10(vMetric)));
    
    % Get the maximum change in the metric
    [max_Delta, idxMaxChange] = max(vDelta);
    
    
end