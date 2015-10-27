% Used in APF
function H1 = BuildH1(m1,m2)
% Build H1, the diagonal matrix of binomials coefficients.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                               Inputs.


% m1 :  Degree of polynomial f with respect to x

% m2 :  Degree of polynomial f with respect to y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise counter
count =1 ;

% get number of diagonals in the matrix
for tot = 0:1:m1+m2+1
    
    for i1 = tot:-1:0
        i2 = tot - i1;
        % if i1 is within the number of rows
        % if i2 is within the number of columns
        if i2 <= m2 && i1 <= m1
            
            % Get binomial coefficients.
            temp_vec(count) = nchoosek(m1,i1) * nchoosek(m2,i2);
            % Increment counter
            count = count+1;
        end
    end
end

H1 = diag(1./temp_vec);

end