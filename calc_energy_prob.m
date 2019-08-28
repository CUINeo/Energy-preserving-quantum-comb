Struct = load('Txe');
T = Struct.T;

C = C_generator(T);
for e = -3 : 3
    p = trace(C * kron(kron(Ux_matrix_gen(), eye(4)), energy_mat_gen(e))) / 4;
    fprintf("Probability of energy %d is %f.\n", e, p);
end

function C = C_generator(T)
    % This function returns the original choi operator
    C = zeros(1792, 1792);
    for x = 1 : 4
        for e = 1 : 7
            C = C + kron(kron(squeeze(T(x, e, :, :)), x_mat_gen(x)), energy_mat_gen(e-4));
        end
    end
end

function matrix = x_mat_gen(x)
    % This function returns the matrix of 4 dimensions with (x, x) = 1
    matrix = zeros(4, 4);
    matrix(x, x) = 1;
end

function matrix = energy_mat_gen(e)
    % This function returns the matrix of 7 dimensions with (e+4, e+4) = 1
    matrix = zeros(7, 7);
    matrix(e+4, e+4) = 1;
end

function Ux_matrix = Ux_matrix_gen()
    % This function returns summation of Ux matrix
    Ux_matrix = zeros(64, 64);
    for x = 0 : 3
        Ux = Ux_generator(x);
        Ux_matrix = Ux_matrix + Ux * ctranspose(Ux);
    end
end

function Ux = Ux_generator(x)
    % This function returns Ux as a ket
    d = 64;
    temp = zeros(8, 8);
    
    for m = 0 : 3
        for n = 0 : 1
            ketM = zeros(4, 1);
            ketM(m+1) = 1;
            
            braN = zeros(1, 2);
            braN(n+1) = 1;
            
            ketN = zeros(2, 1);
            n_ = mod(fx(x, m)+n, 2);
            ketN(n_+1) = 1;
            
            temp = temp + kron(ketM * ctranspose(ketM), ketN * braN);
        end
    end
    
    Ux = reshape(temp, [d 1]);
end

function y = fx(x, m)
    % This function returns the value of fx(m)
    if m == x
        y = 1;
    else
        y = 0;
    end
end