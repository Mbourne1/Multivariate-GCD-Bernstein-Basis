function [fx_elevated] = DegreeElevate_Univariate(fx, r)
% Function performs degree elevation on polynomial f(x).
%
% % Inputs
% 
% fx : (Vector) Coefficients of polynomial f(x)
%
% r : (Int) Number of degree elevations
%
%
% % Outputs
%
% fx_elevated : (Vector) Coefficients of degree elevated f(x)




% Get the degree of polynomial f(x)
m = GetDegree_Univariate(fx);

% Initialise vector of coefficients of degree elevated form
fx_elevated = zeros(m+r+1, 1);

% for each coefficient in the degree elevated polynomial
for k = 0 : 1 : m+r
    
   temp_sum = 0;
   
   for j = max(0,k-r):1:min(m,k)
       temp_sum = temp_sum + (fx(j+1) * nchoosek(m,j) * nchoosek(r,k-j)) ./ nchoosek(m+r,k) ;
   end
   
   fx_elevated(k+1) = temp_sum;
   
end


% Get control points of the input polynomial over unit interval
fx = zeros(m + 1, 1);

for i = 0 : 1 : m
    fx(i + 1) = i ./ m;
end
CP_f = [fx'; fx'];

% Get control points of degree elevated polynomial
fx = zeros(m+r+1,1);

for i = 0:1:(m+r)
    fx(i+1) = i./(m+r);
end

CP_f_elv = [fx'; fx_elevated' ];

CP{1} = CP_f;
CP{2} = CP_f_elv;




end