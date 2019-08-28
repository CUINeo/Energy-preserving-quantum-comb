% version 3 of the full oracle case, with ancilla energy {-3, -2, -1, 0, 1, 2, 3}
% outcome is 0.8125
d1 = 8;
d2 = 8;
n = d1 * d2;
x_num = 4;
e_num = 7;

Omega0 = Omega_generator(0);
Omega1 = Omega_generator(1);
Omega2 = Omega_generator(2);
Omega3 = Omega_generator(3);

Projectorm3 = Projector_generator(-3);
Projectorm2 = Projector_generator(-2);
Projectorm1 = Projector_generator(-1);
Projector0 = Projector_generator(0);
Projector1 = Projector_generator(1);
Projector2 = Projector_generator(2);
Projector3 = Projector_generator(3);

cvx_clear
cvx_begin sdp
cvx_solver sedumi
    variable T(x_num, e_num, n, n)

    S = zeros(n, n);
    for i = 1 : 4
        for j = 1 : 7
            S = S + squeeze(T(i, j, :, :));
        end
    end

    S0 = zeros(n, n);
    S1 = zeros(n, n);
    S2 = zeros(n, n);
    S3 = zeros(n, n);
    for i = 1 : 7
        S0 = S0 + squeeze(T(1, i, :, :));
        S1 = S1 + squeeze(T(2, i, :, :));
        S2 = S2 + squeeze(T(3, i, :, :));
        S3 = S3 + squeeze(T(4, i, :, :));
    end
    maximize(trace(S0 * Omega0) + trace(S1 * Omega1) + trace(S2 * Omega2) + trace(S3 * Omega3))
    
    subject to
        % Comb condition
        trace(S) == d2;
        S == kron(PartialTrace(S, 2, [d1 d2]), eye(d2)/d2);
        for i = 1 : 4
            for j = 1 : 7
                squeeze(T(i, j, :, :)) == hermitian_semidefinite(n);
            end
        end
        % Energy preserving
        for i = 1 : 4
            Projectorm3 * squeeze(T(i, 1, :, :)) * Projectorm3 == squeeze(T(i, 1, :, :));
            Projectorm2 * squeeze(T(i, 2, :, :)) * Projectorm2 == squeeze(T(i, 2, :, :));
            Projectorm1 * squeeze(T(i, 3, :, :)) * Projectorm1 == squeeze(T(i, 3, :, :));
            Projector0 * squeeze(T(i, 4, :, :)) * Projector0 == squeeze(T(i, 4, :, :));
            Projector1 * squeeze(T(i, 5, :, :)) * Projector1 == squeeze(T(i, 5, :, :));
            Projector2 * squeeze(T(i, 6, :, :)) * Projector2 == squeeze(T(i, 6, :, :));
            Projector3 * squeeze(T(i, 7, :, :)) * Projector3 == squeeze(T(i, 7, :, :));
        end
cvx_end

function Omega = Omega_generator(x)
    % This function returns Omega
    d = 64;
    Omega = zeros(d, d);

    Ux = Ux_generator(x);
    Omega = Omega + 1/4 * Ux * ctranspose(Ux);
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

function Projector = Projector_generator(E)
    % This function returns the whole projector with energy difference
    % E1 + E = E2
    n = 64;  
    Projector = zeros(n, n);

    for E1 = 0 : 3
        E2 = E1 + E;
        if E2 >= 0 && E2 <= 3
            P1 = Single_projector_generator(E1);
            P2 = Single_projector_generator(E2);
            Projector = Projector + kron(P1, P2);
        end
    end
end

function SingleProjector = Single_projector_generator(e)
    % This function returns a projector with energy e and dimension 8
    d = 8;
    SingleProjector = zeros(d, d);

    if e == 0
        SingleProjector(1, 1) = 1;
    elseif e == 1
        SingleProjector(2, 2) = 1;
        SingleProjector(3, 3) = 1;
        SingleProjector(5, 5) = 1;
    elseif e == 2
        SingleProjector(4, 4) = 1;
        SingleProjector(6, 6) = 1;
        SingleProjector(7, 7) = 1;
    elseif e == 3
        SingleProjector(8, 8) = 1;
    else
        SingleProjector = eye(d);
    end
end