% reduced oracle case, has some logical errors
d1 = 4;
d2 = 4;
d3 = 4;
n = 64;

Omega = Omega_generator();
Projector = Projector_generator();

cvx_begin sdp
    variable S(n, n) complex semidefinite
    maximize(trace(S*Omega))
    subject to
        % Comb condition
        PartialTrace(S, 3, [d1 d2 d3]) == kron(PartialTrace(S, [2 3], [d1 d2 d3]), eye(d2) / d2);
        trace(S) == d2;
        % Energy preserving
        Projector * S * Projector == S;
cvx_end

function Omega = Omega_generator()
    % This function returns Omega
    n = 64;
    Omega = zeros(n, n);
    
    for x = 1 : 4
        vecX = zeros(4, 1);
        vecX(x) = 1;
        Vx = eye(4);
        Vx(x, x) = -1;
        Omega = Omega + 1/4 * kron(transpose(reshape(Vx, [16 1]) * transpose(reshape(Vx, [16 1]))), vecX * transpose(vecX));
    end
end

function Projector = Projector_generator()
    % This function returns the whole projector
    Projector = zeros(64, 64);
    
    % Total energy is 0
    P2 = Single_projector_generator(0);
    % 00
    P1 = Single_projector_generator(0);
    P3 = Single_projector_generator(0);
    Projector = Projector + kron(kron(P1, P3), P2);
    
    % Total energy is 1
    P2 = Single_projector_generator(1);
    % 10
    P1 = Single_projector_generator(1);
    P3 = Single_projector_generator(0);
    Projector = Projector + kron(kron(P1, P3), P2);
    % 01
    P1 = Single_projector_generator(0);
    P3 = Single_projector_generator(1);
    Projector = Projector + kron(kron(P1, P3), P2);
    
    % Total energy is 2
    P2 = Single_projector_generator(2);
    % 11
    P1 = Single_projector_generator(1);
    P3 = Single_projector_generator(1);
    Projector = Projector + kron(kron(P1, P3), P2);
    % 20
    P1 = Single_projector_generator(2);
    P3 = Single_projector_generator(0);
    Projector = Projector + kron(kron(P1, P3), P2);
    % 02
    P1 = Single_projector_generator(0);
    P3 = Single_projector_generator(2);
    Projector = Projector + kron(kron(P1, P3), P2);
    
    % Total energy is 3
    P2 = Single_projector_generator(3);
    % 03
    P1 = Single_projector_generator(0);
    P3 = Single_projector_generator(3);
    Projector = Projector + kron(kron(P1, P3), P2);
    % 12
    P1 = Single_projector_generator(1);
    P3 = Single_projector_generator(2);
    Projector = Projector + kron(kron(P1, P3), P2);
    % 21
    P1 = Single_projector_generator(2);
    P3 = Single_projector_generator(1);
    Projector = Projector + kron(kron(P1, P3), P2);
    % 30
    P1 = Single_projector_generator(3);
    P3 = Single_projector_generator(0);
    Projector = Projector + kron(kron(P1, P3), P2);
end

function SingleProjector = Single_projector_generator(e)
    % This function returns a projector with energy e
    d = 4;
    
    SingleProjector = zeros(d, d);
    SingleProjector(e+1, e+1) = 1;
end