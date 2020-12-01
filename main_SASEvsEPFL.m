% Script to generate results on EPFL campus 5 node distribution grid

close all
clear
clc

addpath(genpath('./'));

% load DATA structure with data from EPFL campus
load('./data/dataEPFL_2014_Nov_17-10-11h.mat')

% set simulation parameters
PMURate = 50; % [Hz]
FromMin = 0;
FromSec = 200;
SimTime = 6; % [sec]
SelectedPhase = 1; % 
LabelPMU = [1 2 3 4 5]; % label of nodes where PMU are located
LabelPMUMeas = [1 2 3]; % label of nodes where PMU are used for estimation
LabelPQ = [2 3 4 5 6];
Label0Inj = 6; % label of zero injection node in the distribution grid
LabelNodes = [1 2 3 4 5 6];
NumPMU = length(LabelPMU); % number of PMUs

sigma_off  = 2e-4; % [rad] standard deviation of sync error offset parameter (beta)
sigma_skew = 1e-2; % [rad] standard deviation of sync error skew parameter (alpha)

% load grid case and base scaling (non in p.u.)
epfl5_grid

% support quantities
MeasIdx = (FromMin * 60 + FromSec) * PMURate + (1:PMURate * SimTime)';
Delta = repmat(1/PMURate * (1:PMURate)', SimTime, 1);

% extract measurement of interest in Volt, Ampere and Watt
[V, I, S, DeltaPhase, SE] = subsetEPFLdata(NumPMU, MeasIdx, SelectedPhase, G.PCC, DATA);

% normalize quantities
V = V./G.baseV;
Vpcc = mean(abs(V(:,G.PCC)));
V = V./Vpcc;
I = I./G.baseA;
S = S./G.baseVA./Vpcc;
G.L = G.L*G.baseV*Vpcc/G.baseA;
G.X = G.X/G.baseV/Vpcc*G.baseA;
G.U0 = G.U0/G.baseV;


%% estimation
% Kalman Filter based algorithm according to paper 
%   M. Pignati, et al. Real-Time State Estimation of the EPFL-Campus Medium-Voltage Grid by Using PMUs. 
%   The Sixth Conference on Innovative Smart Grid Technologies (ISGT2015), 2015.

% building necessary matrices and running filter
[F, H, R, y, x0, P0] = build_KF_EPFL_matrices(LabelPMUMeas, Label0Inj, G, V, I, DeltaPhase);
[vEPFL, thetaEPFL, Pepfl] = KalmanFilterEPFL(F, H, R, y, x0, P0);

% Proposed Kalman Filter  (SASE) algorithm
measurement_type = 'voltage'; % 'voltage+power' % define which type of measurements to use
[F, Q, H, R, y, x0, P0, u0, s0] = build_KF_SASE_matrices(...
    LabelPMUMeas, Label0Inj, G, V, I, S, sigma_off, sigma_skew, Delta, measurement_type);

[xSASE, Psase] = KalmanFilterSASE(F, Q, H, R, y, x0, P0);
[vSASE, thetaSASE, vSASE_pf, thetaSASE_pf, offSASE, skewSASE, ...
    PvSASE, PthetaSASE, PoffSASE, PskewSASE, pSASE, qSASE, PpSASE, ...
        PqSASE] = extract_KF_SASE_info(LabelPMUMeas, Label0Inj, G, V, I, S, xSASE, Psase);


%% postprocess results    
% extract phase PCC
VApcc = unwrap(angle(V(:, G.PCC)));

% angles w.r.t. PCC
dthetaEPFL = thetaEPFL - repmat(VApcc, [1, size(thetaEPFL, 2)]);

% additional electric quantities
uEPFL = vEPFL .* exp(1i*dthetaEPFL);
iEPFL = uEPFL * G.L;
sEPFL = uEPFL .* conj(iEPFL);
pEPFL = real(sEPFL);
qEPFL = imag(sEPFL);

sSASE = pSASE + 1i.*qSASE;

% TVEs
% do not consider the transient for the EPFL so that it is at steady state
Idx = PMURate + (1:(SimTime - 1) * PMURate);

TVE_s_epfl = computeTVE(sEPFL(Idx, LabelPMUMeas), S(Idx, LabelPMUMeas));
TVE_s_sase = computeTVE(sSASE(Idx, LabelPMUMeas), S(Idx, LabelPMUMeas));
TVE_s_pred_epfl = computeTVE(sEPFL(Idx, setdiff(LabelPMU, LabelPMUMeas)), S(Idx, setdiff(LabelPMU, LabelPMUMeas)));
TVE_s_pred_sase = computeTVE(sSASE(Idx, setdiff(LabelPMU, LabelPMUMeas)), S(Idx, setdiff(LabelPMU, LabelPMUMeas)));


%% print and plot results
fprintf(1,'\n TVE Power EPFL:\t\t %d', TVE_s_epfl);
fprintf(1,'\n TVE Power SASE:\t\t %d', TVE_s_sase);
fprintf(1,'\n TVE Power EPFL (prediction):\t %d', TVE_s_pred_epfl);
fprintf(1,'\n TVE Power SASE (prediction):\t %d', TVE_s_pred_sase);
fprintf(1,'\n');

% uncomment to run
%plot_SASEvsEPFL

