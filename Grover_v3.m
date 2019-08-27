% version 3 of the full oracle case, inaccurate output
d1 = 8;
d2 = 8;
n = d1 * d2;

Omega0 = Omega_generator(0);
Omega1 = Omega_generator(1);
Omega2 = Omega_generator(2);
Omega3 = Omega_generator(3);

Projector0 = Projector_generator(0);
Projector1 = Projector_generator(1);
Projector2 = Projector_generator(2);
Projector3 = Projector_generator(3);

cvx_clear
cvx_begin sdp
cvx_solver sedumi
    variable T00(n, n) complex semidefinite
    variable T01(n, n) complex semidefinite
    variable T02(n, n) complex semidefinite
    variable T03(n, n) complex semidefinite
    variable T10(n, n) complex semidefinite
    variable T11(n, n) complex semidefinite
    variable T12(n, n) complex semidefinite
    variable T13(n, n) complex semidefinite
    variable T20(n, n) complex semidefinite
    variable T21(n, n) complex semidefinite
    variable T22(n, n) complex semidefinite
    variable T23(n, n) complex semidefinite
    variable T30(n, n) complex semidefinite
    variable T31(n, n) complex semidefinite
    variable T32(n, n) complex semidefinite
    variable T33(n, n) complex semidefinite
    
    S0 = T00 + T01 + T02 + T03;
    S1 = T10 + T11 + T12 + T13;
    S2 = T20 + T21 + T22 + T23;
    S3 = T30 + T31 + T32 + T33;
    maximize(trace(S0 * Omega0) + trace(S1 * Omega1) + trace(S2 * Omega2) + trace(S3 * Omega3))
    
    subject to
        % Comb condition
        trace(T00+T01+T02+T03+T10+T11+T12+T13+T20+T21+T22+T23+T30+T31+T32+T33) == d2;
        T00+T01+T02+T03+T10+T11+T12+T13+T20+T21+T22+T23+T30+T31+T32+T33 == kron(PartialTrace(T00+T01+T02+T03+T10+T11+T12+T13+T20+T21+T22+T23+T30+T31+T32+T33, 2, [d1 d2]), eye(d2) / d2);
        % Energy preserving
        Projector0 * T00 * Projector0 == T00;
        Projector0 * T10 * Projector0 == T10;
        Projector0 * T20 * Projector0 == T20;
        Projector0 * T30 * Projector0 == T30;
        Projector1 * T01 * Projector1 == T01;
        Projector1 * T11 * Projector1 == T11;
        Projector1 * T21 * Projector1 == T21;
        Projector1 * T31 * Projector1 == T31;
        Projector2 * T02 * Projector2 == T02;
        Projector2 * T12 * Projector2 == T12;
        Projector2 * T22 * Projector2 == T22;
        Projector2 * T32 * Projector2 == T32;
        Projector3 * T03 * Projector3 == T03;
        Projector3 * T13 * Projector3 == T13;
        Projector3 * T23 * Projector3 == T23;
        Projector3 * T33 * Projector3 == T33;
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
    d1 = 8;
    d2 = 8;
    n = 64;
    
    Projector = zeros(n, n);
    
    % E1 = 0
    P1 = Single_projector_generator(0, d1);
    P2 = Single_projector_generator(E, d2);
    Projector = Projector + kron(P1, P2);
    
    % E1 = 1
    if E+1 <= 3
        P1 = Single_projector_generator(1, d1);
        P2 = Single_projector_generator(E+1, d2);
        Projector = Projector + kron(P1, P2);
    end
    
    % E1 = 2
    if E+2 <= 3
        P1 = Single_projector_generator(2, d1);
        P2 = Single_projector_generator(E+2, d2);
        Projector = Projector + kron(P1, P2);
    end
    
    % E1 = 3
    if E+3 <= 3
        P1 = Single_projector_generator(3, d1);
        P2 = Single_projector_generator(E+3, d2);
        Projector = Projector + kron(P1, P2);
    end
end

function SingleProjector = Single_projector_generator(e, d)
    % This function returns a projector with energy e and dimension d
    if d == 2
        SingleProjector = zeros(d, d);
        SingleProjector(e+1, e+1) = 1;
    elseif d == 8
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
    else
        SingleProjector = eye(d);
    end
end