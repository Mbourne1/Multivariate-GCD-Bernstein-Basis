function [implic] = Implicitize_Bezier_Bezout(CP)
% Given a set of Bezier control points for the Bezier curve. Implicitize 
% the Curve. This code corresponds to the implicitization notes in
% Sederberg CAGD notes Section 16.6.
%
% % Inputs.
%
% CP : Vector of control points stored as follows
%        ___________________
%       |x0  x1  x2 ...  xn |
%       |y0  y1  y2 ...  yn |
%       |___________________|
%
% % Outputs.
%
% implic : The coefficients of the polynomial in implicit form.

% Get the degree of the polynomial
n = size(CP,2)-1;

% initialise the matrix to store L_(i,j)
mat = zeros(n,n);
mat = sym(mat);

% For each row
for i=0:1:(n-1)
    % for each column
    for j= 0:1:(n-1)
        
        mat(i+1,j+1) = Build_Lij((n-1)-i,(n-1)-j,CP);
        
    end
end

mat
mydet(mat)

implic = expand(det(mat));


end

function Lij = Build_Lij(i,j,CP)

% Get the degree of the polynomial
n = size(CP,2)-1;

sum = 0;

for m = 0:1:min(i,j)
    k = i+j+1-m;
    
    if k <= n
        % set i = k
        
        
        xi = CP(1,k+1);
        yi = CP(2,k+1);
        xj = CP(1,m+1);
        yj = CP(2,m+1);
        
        
        sum = sum + Build_lij(n,k,m,xi,yi,xj,yj);
    end
end

Lij = sum;

end


function lij = Build_lij(n,i,j,xi,yi,xj,yj)

wi = 1;
wj = 1;


x = sym('x');
y = sym('y');


mat = ...
    [
    x   y   1;
    xi  yi  1;
    xj  yj  1;
    ];

lij = ...
    nchoosek(n,i) * ...
    nchoosek(n,j) * ...
    wi ...
    *   wj ...
    * det(mat);

end