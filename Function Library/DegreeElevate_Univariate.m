function [c] = DegreeElevate_Univariate(ft,r)
% Function performs degree elevation on polynomial f(x).
% Let fx be the coefficients of polynomial f(x) of degree m. 
% r = Number of degree elevations such that the output polynomial is of
% degree m + r.

m = length(ft) -1;

% for each c(i) in the vector
c = zeros(m+r+1,1);

for k = 0:1:m+r
   temp_sum = 0;
   for j = max(0,k-r):1:min(m,k)
       temp_sum = temp_sum + (ft(j+1) * nchoosek(m,j) * nchoosek(r,k-j)) ./ nchoosek(m+r,k) ;
   end
   c(k+1) = temp_sum;
end


% Get control points of the input polynomial over unit interval
fx = zeros(m+1,1);
for i = 0:1:m
    fx(i+1) = i./m;
end
CP_f = [fx'; ft'];

% Get control points of degree elevated polynomial
fx = zeros(m+r+1,1);

for i = 0:1:(m+r)
    fx(i+1) = i./(m+r);
end

CP_f_elv = [fx'; c' ];

CP{1} = CP_f;
CP{2} = CP_f_elv;




end