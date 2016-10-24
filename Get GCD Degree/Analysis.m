
function [max_Delta_MaxMin_RowSum_R,indexMaxChange_RowNorm] = Analysis(vRatio_MaxMin_RowNorm_R)

    vRatio_MaxMin_RowNorm_R = sanitize(vRatio_MaxMin_RowNorm_R);
    % % Analyse Max:Min Row Norms for each subresultant
    % Get the change in the ratios from one subresultant to the next.
    vDelta_MaxMin_RowNorm_R = abs(diff(log10(vRatio_MaxMin_RowNorm_R)));
    % Get the maximum change in rowsum ratio and its index
    [max_Delta_MaxMin_RowSum_R,indexMaxChange_RowNorm] = max(vDelta_MaxMin_RowNorm_R);
end