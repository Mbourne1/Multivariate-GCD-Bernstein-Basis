function [fxy, arr_fxy, arr_hxy, arr_wxy,  m] = Examples_Roots_FromCoefficients_Bivariate(ex_num)
%
% % Inputs
% 
% ex_num : (String) Example number as a string
%
% % Outputs
%
% fxy : (Matrix) Coefficients of f(x,y)
%
% m : (Int) Total degree of f(x,y)

syms x y ;

% Add examples folder to current path
addpath(genpath('../Examples'));

% Get the factor and multiplicity array for polynomial f(x,y)
fxy_root_mult_arr = Roots_Examples_Bivariate(ex_num);

% return the coefficient matrix of f(x,y)
fxy = GetCoefficientsFromSymbolicRoots_Bivariate(fxy_root_mult_arr);


arr_fxy = GetArray_fxy_Exact(fxy_root_mult_arr);
arr_hxy = GetArray_hxy_Exact(fxy_root_mult_arr);
arr_wxy = GetArray_wxy_Exact(fxy_root_mult_arr);

% Get the symbolic polynomial in power basis
symbolic_fxy = GetSymbolicPolyFromSymbolicRoots_Bivariate(fxy_root_mult_arr);

% Get the total degree of f(x,y)
m = double(feval(symengine, 'degree', symbolic_fxy,x));

fprintf('%s : Polynomial f(x,y). \n',mfilename);
display(symbolic_fxy);


end


function [arr_fxy] = GetArray_fxy_Exact(fxy_factor_multiplicity_matrix)
%
%
% fxy_factor_multiplicity_matrix : Matrix containing the symbolic factors of
% f(x) and the multiplicity of the factors.
%
% % Outputs
%
% arr_fx : (Array of Vectors) Each vector contains the coefficients of a
% polynomial f_{i}(x) in the sequence generated in Gauss method of
% factorisation.
%
%

% Get symbolic factors
vFactors = (fxy_factor_multiplicity_matrix(:,1));

% Get multiplicity of factors
vMultiplicity_f0 = double(fxy_factor_multiplicity_matrix(:,2));

% Get max multiplicty
maxMultiplicity = max(vMultiplicity_f0);

% Set number of polynomials in array
nPolys_arr_fx = maxMultiplicity + 1;

% Initialise array
arrMultiplicities = cell(nPolys_arr_fx, 1);
arr_fxy = cell(nPolys_arr_fx, 1);

for i = 1 : 1 : nPolys_arr_fx
    
    % Get multiplicity of factors in f_{i}(x)
    vMultiplicity_fi = vMultiplicity_f0 - (i-1);
    vMultiplicity_fi(vMultiplicity_fi < 0) = 0;
    
    arrMultiplicities{i} = vMultiplicity_fi;
    
    % Get coefficients of f_{i}(x)
    if sum(vMultiplicity_fi) > 0
        arr_fxy{i} = GetCoefficientsFromSymbolicRoots_Bivariate([vFactors vMultiplicity_fi]);
    else
        arr_fxy{i} = 1;
    end
end


end




function [arr_hxy] = GetArray_hxy_Exact(fxy_factor_multiplicity_matrix)
%
% % Inputs
%
% fx_factor_multiplicity_matrix : (Matrix) containing symbolic factors and
% their multiplicities
%


% Get multiplicity structure of f_{0}(x)
vMultiplicity_f0 = double(fxy_factor_multiplicity_matrix(:,2));

% Get the set of factors of f_{0}(x)
vFactors = (fxy_factor_multiplicity_matrix(:,1));

% Get multiplicity structure of f_{i}(x)
arr_vMultiplicity_fx = GetMultiplicityArr_fx(vMultiplicity_f0);

nPolysArr_fx = max(vMultiplicity_f0) + 1;

% Get number of polynomials h_{i}(x)
nPolys_arr_hx = nPolysArr_fx - 1;

% Initialise array to store h_{i}
arr_hxy = cell(nPolys_arr_hx, 1);

for i = 1 : 1 : nPolys_arr_hx
   
    vMultiplicity_hi = arr_vMultiplicity_fx{i} - arr_vMultiplicity_fx{i+1};
    
    
    arr_hxy{i} = GetCoefficientsFromSymbolicRoots_Bivariate([vFactors vMultiplicity_hi]);
    
    
end



end

function [arr_wxy] = GetArray_wxy_Exact(fxy_factor_multiplicity_matrix)
%
% % Inputs
%
% fx_factor_multiplicity_matrix : (Matrix) containing symbolic factors and
% their multiplicities
%


% Get multiplicity structure of f_{0}(x)
vMultiplicity_f0 = double(fxy_factor_multiplicity_matrix(:,2));

% Get the set of factors of f_{0}(x)
vFactors = (fxy_factor_multiplicity_matrix(:,1));

% Get multiplicity structure of f_{i}(x)
arr_vMultiplicity_fx = GetMultiplicityArr_fx(vMultiplicity_f0);

nPolysArr_fx = max(vMultiplicity_f0) + 1;

% Get number of polynomials h_{i}(x)
nPolys_arr_hx = nPolysArr_fx - 1;

% Get number of polynomials w_{i}(x)
nPolys_arr_wx = nPolys_arr_hx;


% Initialise array to store h_{i}
arr_wxy = cell(nPolys_arr_hx, 1);



% Get maximum multiplicity
maxMultiplicity = max(vMultiplicity_f0);

% Initialise a cell array, each cell contains a vector of symbolic factors.
arr_vFactors = cell(maxMultiplicity,1);

nFactors = size(fxy_factor_multiplicity_matrix, 1);

for i = 1 : 1 : nFactors
   
    % Get factor
    myFactor = fxy_factor_multiplicity_matrix(i,1);
    
    % Get multiplicity of factor
    myMultiplicity = fxy_factor_multiplicity_matrix(i,2);
    
    % Get vector of factors of multiplicity 'myMultiplicity' 
    vFactors = arr_vFactors{myMultiplicity};
    vFactors = [vFactors ; myFactor];
    arr_vFactors{myMultiplicity} = vFactors;
    
end



for i = 1 : 1 : maxMultiplicity
   
    vFactors = arr_vFactors{i};
    vMultiplicity = ones(length(vFactors),1);
    
    if length(vFactors) >= 1
        arr_wxy{i} = GetCoefficientsFromSymbolicRoots_Bivariate([vFactors vMultiplicity]);
    else
        arr_wxy{i} = 1;
    end
end






end


function [arrMultiplicities] = GetMultiplicityArr_fx(vMultiplicity_f0)
%
%
% fx_factor_multiplicity_matrix : Matrix containing the symbolic factors of
% f(x) and the multiplicity of the factors.
%
% % Outputs
%
% arr_fx : (Array of Vectors) Each vector contains the coefficients of a
% polynomial f_{i}(x) in the sequence generated in Gauss method of
% factorisation.


% Get max multiplicty
maxMultiplicity = max(vMultiplicity_f0);

% Set number of polynomials in array
nPolys_arr_fx = maxMultiplicity + 1;

% Initialise array
arrMultiplicities = cell(nPolys_arr_fx, 1);


for i = 1 : 1 : nPolys_arr_fx
    
    % Get multiplicity of factors in f_{i}(x)
    vMultiplicity_fi = vMultiplicity_f0 - (i-1);
    vMultiplicity_fi(vMultiplicity_fi < 0) = 0;
    
    arrMultiplicities{i} = vMultiplicity_fi;
    
    
    
end


end