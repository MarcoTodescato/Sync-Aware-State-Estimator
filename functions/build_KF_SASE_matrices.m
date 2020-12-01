function [F, Q, H, R, y, x0, P0, u0, s0] = build_KF_SASE_matrices(...
    LblPMU, Lbl0Inj, G, V, I, S, sigma_off, sigma_skew, Delta, meas_type)
% Computes all the quantities necessary to run the SASE filter where the 
% state vector is expressedn in powers. No reference node is chosen.
%
% Args:
%   LblPMU: label of nodes where PMU measurements are take from
%   Lbl0Inj: label of zero-injection node
%   G: struct with grid
%   V: complex voltages measurements
%   I: complex currents measurements
%   S: complex powers measurements
%   sigma_off: std of sync offset error parameter (beta)
%   sigma_skew: std of sync skew error parameter (alpha)
%   Delta: vector of discrete time instant PMU measurements are taken
%       during two consequitve GPS sync signals
%   meas_type: type of measurements used
%
% Returns:
%   F, Q, H, R: matrix for Kalman Filter
%   y: measurements
%   x0, P0: initial condition for state and covariance
%   u0: nominal working point for voltage
%   s0: nominal working point for power

PCC = G.PCC;
L = G.L;
baseV = G.baseV;
baseA = G.baseA;

TotNumNodes = size(L,1);
[NumIter, NumNodes] = size(V);

LblNo0Inj = setdiff(1:TotNumNodes, Lbl0Inj);
Lno0inj = L(LblNo0Inj, LblNo0Inj) - L(LblNo0Inj, Lbl0Inj) * ( L(Lbl0Inj, Lbl0Inj) \ L(Lbl0Inj, LblNo0Inj));

idx_pmu = logical(pad(NumNodes,LblPMU,1));
idx_pmu2 = logical(kron(idx_pmu,ones(2,1)));

NumPMU = length(LblPMU);
StateDim = 2*(NumNodes + NumPMU);

% normalize measurements wrt PCC
Vpcc = mean(abs(V(:,PCC))); % pcc mean absolute value
V = V ./ Vpcc;
S = S ./ Vpcc;    % scale powers
Lno0inj = Lno0inj .* Vpcc;

% compute average/predicited/nominal working point
s0 = [mean(S,1).' ; 0];
[u0, i0] = networkState(G, s0, zeros(TotNumNodes, 1));

% linearization matrix
Id = eye(NumNodes);
Hdelay2v = c2ri(1i*eye(NumPMU),'evol');
Hs2s = c2ri(Id(idx_pmu,:));
Ru = ri2p_jacobian(u0(LblNo0Inj));
N = kron(Id,[1 0 ; 0 -1]);
Hv2s = (c2ri(diag(conj(i0(LblNo0Inj)))) + c2ri(diag(u0(LblNo0Inj)))*N*c2ri(Lno0inj))*Ru;
Hs2v = inv(Hv2s);
Hs2vPMU = Hs2v(idx_pmu2,:);

% initial conditions
x0 = zeros(StateDim, 1);
P0 = blkdiag(1e-11 * eye(2*NumNodes),...
    sigma_off^2 * eye(NumPMU),...
    sigma_skew^2 * eye(NumPMU));

% build state and output matrices
F = zeros(StateDim,StateDim,NumIter);
Q = zeros(StateDim,StateDim,NumIter);
for t = 1:NumIter
    
    % build state matrix
    if t > 1 && Delta(t) == Delta(1)
        F(:,:,t) = blkdiag(eye(2*NumNodes),...
                    [eye(NumPMU)   , Delta(t-1)*eye(NumPMU); ...
                     zeros(NumPMU) , eye(NumPMU)]);
    else
        F(:,:,t) = eye(StateDim);
    end
    
    % build process covariance matrix
    if t>1 && Delta(t) == Delta(1)
        Q(:,:,t) = P0/100;
    end
end

% build measurements, output matrix and measurements correlation matrix
y_v = zeros(2*NumPMU,NumIter);
y_v(1:2:end,:) = abs(V(:,LblPMU))';
y_v(2:2:end,:) = unwrap(angle(V(:,LblPMU)))';
dy_v = y_v - repmat(c2p(u0(LblPMU)),[1,NumIter]);

R = computeCovMatrixEPFLMeasurementsSASE(LblPMU,PCC,baseV*Vpcc,baseA,V,I);

switch meas_type 
    
    case 'voltage+power'
               
        y_s = zeros(2*NumPMU,NumIter);
        y_s(1:2:end,:) = real(S(:,LblPMU).');
        y_s(2:2:end,:) = imag(S(:,LblPMU).');
        dy_s = y_s - repmat(c2ri(s0(LblPMU),'v'),[1,NumIter]);
        
        y = [dy_v;dy_s];
        
        NumMeas = size(y,1);
                
        H = zeros(NumMeas,StateDim,NumIter);
        for t = 1:NumIter
            % output matrix
            H(:,:,t) = [Hs2vPMU , [Hdelay2v , Hdelay2v*Delta(t)]; ...
                        Hs2s , zeros(2*NumPMU)];
        end
        
    otherwise % only voltages
        
        R = R(1:2*NumPMU,1:2*NumPMU,:);
        
        y = dy_v;
        
        NumMeas = size(y,1);
        
        H = zeros(NumMeas,StateDim,NumIter);
        for t = 1:NumIter
            % output matrix
            H(:,:,t) = [Hs2vPMU , [Hdelay2v , Hdelay2v*Delta(t)]];
        end
end                
end

