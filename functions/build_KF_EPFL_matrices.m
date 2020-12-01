function [F, H, R, y, x0, P0] = build_KF_EPFL_matrices(LblPMU, Lbl0Inj, G, V, I, DeltaPhase)
% Build necessary matrices for the KF-based algorithm proposed in
%   M. Pignati, et al. Real-Time State Estimation of the EPFL-Campus Medium-Voltage Grid by Using PMUs. 
%   The Sixth Conference on Innovative Smart Grid Technologies (ISGT2015), 2015.
%
% Args:
%   LblPMU: label of node where PMUs measurements are used for estimation
%   Lbl0Inj: label of zero-injection node
%   G: grid struct
%   V: voltage measurements
%   I: current measurements
%   DeltaPhase: cumulative phase measurements
%
% Returns:
%   F, H, R matrices for the linear model used by the KF
%   y: measurements
%   x0: initial condition
%   P0: initial covariance

PCC = G.PCC;
L = G.L;
baseV = G.baseV;
baseA = G.baseA;
                                     
NumNodes = size(L, 1);
NumIter = size(V, 1);

% useful logical indexing vectors
idx_pmu = logical(pad(NumNodes, LblPMU, 1));
idx_0inj = logical(pad(NumNodes, Lbl0Inj, 1));
idx_pmu2 = logical(kron(idx_pmu, [1; 1]));
idx_all = logical(kron(idx_pmu | idx_0inj, [1; 1]));

% build state matrix
Id = eye(NumNodes);
theta = [0 ; diff(DeltaPhase)];
F = zeros(2*NumNodes,2*NumNodes,NumIter);
for t = 1:NumIter
    
    theta_t = theta(t,:); 
    
    F(:,:,t) = [cos(theta_t) * Id , -sin(theta_t) * Id; ...
                sin(theta_t) * Id ,  cos(theta_t) * Id];
        
end

% build output matrix
H = [c2ri(Id(idx_pmu, :));
     c2ri(L(idx_pmu | idx_0inj, :))];

% build measurements in real-imag (consider zero injection if necessary)
y_v  = zeros(2 * NumNodes, NumIter);
y_v(2 * LblPMU - 1, :) = real(V(:, LblPMU))';
y_v(2 * LblPMU, :) = imag(V(:, LblPMU))';

y_i = zeros(2 * NumNodes, NumIter);
y_i(2 * LblPMU - 1, :) = real(I(:, LblPMU))';
y_i(2 * LblPMU, :) = imag(I(:, LblPMU))';

y = [y_v(idx_pmu2, :); y_i(idx_all, :)];
    
% measurement covariance matrix
NumMeas = size(y, 1);
R = zeros(NumMeas, NumMeas, NumIter);
r_v = zeros(2 * NumNodes, 1);
r_i = zeros(2 * NumNodes, 1);

for t = 1:NumIter
    
    % extract current meas
    v = V(t, LblPMU).';
    i = I(t, LblPMU).';
    
    %[std_vm,std_va,std_im,std_ia,...
    [~, ~, ~, ~, std_v_re, std_v_im, std_i_re, std_i_im] = computeStdMeasEPFL(v, i, LblPMU, PCC, baseV, baseA);
 
    r_v(2 * LblPMU - 1) = std_v_re.^2;
    r_v(2 * LblPMU) = std_v_im.^2;
    r_i(2 * LblPMU - 1) = std_i_re.^2; 
    r_i(2 * LblPMU) = std_i_im.^2;
    
    % zero inj std
    r_i(2 * find(idx_0inj) - 1) = min(std_i_re).^2;
    r_i(2 * find(idx_0inj))   = min(std_i_im).^2;
    
    R(:, :, t) = diag([r_v(idx_pmu2); r_i(idx_all)]); 
end

% initial conditions
x0 = zeros(2 * NumNodes, 1);
P0 = 1e2 * eye(2 * NumNodes); 

end

