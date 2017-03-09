function GM_fx = GetGeometricMean_MyMethod(fxy,n1_k1,n2_k2)
% GetGeometricMean_MyMethod(fxy,n1_k1,n2_k2)
%
% Inputs.
%
% fxy : (Matrix) Coefficients of polynomial f(x,y)
%
% n1_k1 : (Int) Degree of polynomial v(x,y) with respect to x
%
% n2_k2 : (Int) Degree of polynomial v(x,y) with respect to y
%
% Outputs.
%
% GM_fx : Geometric mean of non-zero entries of the first partition of the
% Sylvester subresultant matrix S_{k_{1},k_{2}}


% Get the degree of polynomial f(x,y)
[m1, m2] = GetDegree_Bivariate(fxy);

% Get the geometric mean of the coefficient parts only
prod_a = 1;
sum_log_a = 0;

for i1 = 0:1:m1
    for i2 = 0:1:m2
        sum_log_a = sum_log_a + log10(fxy(i1+1,i2+1));
    end
end


for i1 = 0:1:m1
    for i2 = 0:1:m2
        prod_a = prod_a .* fxy(i1+1,i2+1);
    end
end


geo_mean_a = prod_a .^ (1./((m1+1)*(m2+1)));

geo_mean_log_a = (1./((m1+1)*(m2+1))) * sum_log_a;
geo_mean_a_2 = 10.^(geo_mean_log_a);



% %
% %
% %
% Get geometric mean of the binomial nchoosek(i1+j1,i1)
prod_b = 1;
for i1 = 0:1:m1
    for j1 = 0:1:n1_k1
        prod_b = prod_b .* nchoosek(i1 + j1 , i1);
    end
end

geo_mean_b = prod_b .^ (1./ ((m1+1)*(n1_k1+1)) );

% %
% %
% %
% Get geometric mean of the binomial nchoosek(i2+j2,i2)
prod_c = 1;
for i2 = 0:1:m2
    for j2 = 0:1:n2_k2
        prod_c = prod_c .* nchoosek(i2 + j2 , i2);
    end
end

geo_mean_c = prod_c .^ (1./ ((m2+1)*(n2_k2+1)));

% %
% %
% %
% Compute geometric mean

GM_fx = geo_mean_a * geo_mean_b * geo_mean_b * geo_mean_c * geo_mean_c ...
    ./ (nchoosek(m1+n1_k1,n1_k1) * nchoosek(m2+n2_k2,n2_k2)); 



end