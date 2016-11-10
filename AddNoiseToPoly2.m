function [fxy_noisy, noise_matrix] = AddNoiseToPoly2(fxy_matrix,el,em)

% Initialise global variables
global SETTINGS

% Set seed
rng(SETTINGS.SEED)



switch SETTINGS.BOOL_NOISE
    case 'y'
        % Get the degrees of polynomial f(x,y)
        [m1,m2] = GetDegree(fxy_matrix);
        
        
        % Obtain a matrix of randomly distributed numbers between [-1 and 1]
        rp = (2*rand(m1+1,m2+1))-ones(m1+1,m2+1);
        
        % multiply by the noise:signal ratio so so that we have a range of
        % values between [-noise : + noise]
        s = rp*el;
        
        % Get the noise matrix
        noise_matrix = fxy_matrix.*s;
        
        % add the noise matrix to the exact matrix
        fxy_noisy = fxy_matrix + noise_matrix;
    case 'n'
        fxy_noisy = fxy_matrix;
        noise_matrix = zeros(size(fxy_matrix));
    otherwise
        error('Error: SETTINGS.BOOL_NOISE')
end
end