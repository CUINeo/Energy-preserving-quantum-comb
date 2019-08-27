% version 1 of the full oracle case, cannot run due to lack of memory
d1 = 8;
d2 = 8;
d3 = 4;
d4 = 2;
n = d1 * d2 * d3 * d4;

Omega = Omega_generator();
Projector = Projector_generator();

cvx_clear
cvx_begin sdp
cvx_solver sedumi
    variable S(n, n) complex semidefinite
    maximize(trace(S * Omega))
    subject to
        % Comb condition
        PartialTrace(S, [3 4], [d1 d2 d3 d4]) == kron(PartialTrace(S, [2 3 4], [d1 d2 d3 d4]), eye(d2) / d2);
        trace(S) == d2;
        % Energy preserving
        Projector * S * Projector == S;
cvx_end

function Omega = Omega_generator()
    % This function returns Omega
    d = 512;
    Omega = zeros(d, d);
    
    for x = 1 : 4
        vecX = zeros(4, 1);
        vecX(x) = 1;
        
        Ux = Ux_generator(x);
        Omega = Omega + 1/4 * kron(kron(Ux * ctranspose(Ux), vecX * ctranspose(vecX)), eye(2));
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

function Projector = Projector_generator()
    % This function returns the whole projector
    d1 = 8;
    d2 = 8;
    d3 = 4;
    d4 = 2;
    n = 512;
    
    I3 = eye(d3);
    Projector = zeros(n, n);
    
    % Total energy is 0
    % 00
    P1 = Single_projector_generator(0, d1);
    P2 = Single_projector_generator(0, d2);
    P4 = Single_projector_generator(0, d4);
    Projector = Projector + kron(P1, kron(I3, kron(P4, P2)));
    
    % Total energy is 1
    % 01
    P1 = Single_projector_generator(0, d1);
    P2 = Single_projector_generator(1, d2);
    P4 = Single_projector_generator(1, d4);
    Projector = Projector + kron(P1, kron(I3, kron(P4, P2)));
    % 10
    P1 = Single_projector_generator(1, d1);
    P2 = Single_projector_generator(1, d2);
    P4 = Single_projector_generator(0, d4);
    Projector = Projector + kron(P1, kron(I3, kron(P4, P2)));
    
    % Total energy is 2
    % 20
    P1 = Single_projector_generator(2, d1);
    P2 = Single_projector_generator(2, d2);
    P4 = Single_projector_generator(0, d4);
    Projector = Projector + kron(P1, kron(I3, kron(P4, P2)));
    % 11
    P1 = Single_projector_generator(1, d1);
    P2 = Single_projector_generator(2, d2);
    P4 = Single_projector_generator(1, d4);
    Projector = Projector + kron(P1, kron(I3, kron(P4, P2)));
    
    % Total energy is 3
    % 30
    P1 = Single_projector_generator(3, d1);
    P2 = Single_projector_generator(3, d2);
    P4 = Single_projector_generator(0, d4);
    Projector = Projector + kron(P1, kron(I3, kron(P4, P2)));
    % 21
    P1 = Single_projector_generator(2, d1);
    P2 = Single_projector_generator(3, d2);
    P4 = Single_projector_generator(1, d4);
    Projector = Projector + kron(P1, kron(I3, kron(P4, P2)));
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