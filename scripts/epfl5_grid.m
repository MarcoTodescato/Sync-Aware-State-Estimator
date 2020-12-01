% Struct representation of the EPFL 5 node distribution test-bed

% number of nodes
G.n = 6;

% incidence matrix
A = horzcat( ...
  ee([-1 2], G.n), ...
  ee([-2 3], G.n), ...
  ee([-3 4], G.n), ...
  ee([-4 5], G.n), ...
  ee([-5 6], G.n));

G.A = A';

% number of edges
G.m = size(G.A,1);

LinesLength = [.460, .073, .072, .035, .311]';

Resistance_xKm = .159; % ohm
Reactance_xKm = .113; % ohm
LineChargSusceptance_xKm = 84.8 * 1e-6; % siemens

% vector of impedances
R = LinesLength * Resistance_xKm;
X = LinesLength * Reactance_xKm;
B = LinesLength * LineChargSusceptance_xKm;

G.z = R + 1i*X;
G.b = B;

% index of the PCC
G.PCC = 1;

% nominal voltage
G.U0 = 11547;

% base values
G.baseMVA = 100;
G.baseKV = G.U0/1e3;
G.baseKA = G.baseMVA/G.baseKV;
G.baseZ = G.baseKV/G.baseKA;
G.baseY = 1/G.baseZ;

G.baseVA = G.baseMVA*1e6;
G.baseV = G.baseKV*1e3;
G.baseA = G.baseKA*1e3;

% vector of static model exponents
G.eta = [0, 0, 0, 0, 0, 0]';
    
% weighted laplacian
G.L = G.A' * diag(1./G.z) * G.A + diag(diag(G.A' * diag(1i*G.b/2) * G.A));

% Green-like matrix
X = inv([G.L ee(1,G.n); ee(1,G.n)' 0]);
G.X = X(1:G.n, 1:G.n);

clear A X R;
