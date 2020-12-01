% Main script to test the proposed SASE algorithm on synthetic testcases
% 
% The scipt can be used to generate results on two different testcases
% using two different execution modes:
%   1. as function of the measurement collected by a fixed number of
%   deployed PMUs
%   2. as function of the number of PMU deployed collecting a fixed number
%   of measurement each

close all
clear
clc

addpath(genpath('./'));

% select execution mode
mode_1 = 'performance_as_function_of_number_of_pmu_measurements';
mode_2 = 'performance_as_function_of_number_of_pmu_deployed';

exe_mode = mode_1;

% toogle save data
saveData = 0;

% data parameters from MATPOWER
define_constants;
P.msg_non_converged = 'SOLVER convergence error: loadflow not feasible';
P.mpopt = mpoption('verbose',0,'out.all',0);

%matpower_case = 'case15';   % indian distribution
matpower_case = 'case123';  % ieee 123 
G.mpc = loadcase(matpower_case); 

% parameters for generating prior and measurements data
P.solver_type = 'matpower';
P.meas_model = 'true';
P.coordinate_type = 'polar'; % estimation algorithm works only for polar coordinates
P.state_correlation = 'correlated';

P.tau = 1;    % temporal window [sec] where we take measures
P.NumWindow = 1; % number of windows (time intervals) measurements are taken for
% load positions (node label) where PMUs are deployed
switch matpower_case
    case 'case123'
        load('data/PMU_Positions_IEEE123.mat');
        P.PMU_POS = PMU_POS;
        clear PMU_POS;
        PmuStep = 5;
    otherwise % case15 
        P.PMU_POS = [3, 7, 13, 15, 10, 14, 8, 12, 5, 11, 6, 9, 4, 2]';
        PmuStep = 1;
end

switch exe_mode
    case mode_1
        P.NumMC = 500;
        P.MaxNumMeas = 60;
        P.NUM_MEAS = [1 2 4 8 15 20 25 30 50 60];
        P.MaxNumPmu = length(P.PMU_POS);
        P.MaxNumIter = P.MaxNumMeas * P.NumWindow; % max number of algorithmic iterations
        P.NUM_PMU = 60;
        save_suffix = '_FuncNumPmuMeasurements';
        
    case mode_2
        P.NumMC = 1; %500
        P.MaxNumMeas = 3; %30
        P.NUM_MEAS = P.MaxNumMeas;
        P.MaxNumPmu = length(P.PMU_POS);
        P.MaxNumIter = P.MaxNumMeas*P.NumWindow;
        P.NUM_PMU = 0:PmuStep:P.MaxNumPmu;
        save_suffix = '_FuncNumPmuDeployed';
        
    otherwise
        error('Selected exe_mode not supported')
end

P.sigma_off  = 2e-4; % sync error offset parameter std (beta)
P.sigma_skew = 1e-2; % clock skew error parameter std (alpha)  
P.sigma_pmu  = 1e-3; % std for PMU measurements
P.real_prior_load_variation = 0.5; % fractional variation with respect to 1-day ahead forecast
P.imag_prior_load_variation = 0.5; % fractional variation with respect to 1-day ahead forecast
P.real_percent_load_variation = 0.01;
P.imag_percent_load_variation = 0.01;

% grid parameters
G.n = size(G.mpc.bus,1);
G.PCC = find(G.mpc.bus(:,BUS_TYPE) == 3);
G.idx_pq = ~logical(pad(G.n,G.PCC,1));
G.num_pq = sum(G.idx_pq);
G.v0 = G.mpc.bus(G.PCC,VM);
G.L = full(makeYbus(G.mpc));
X = inv([G.L ee(G.PCC,G.n); ee(G.PCC,G.n)' 0]);
G.X = X(1:G.n,1:G.n);
clear X
G.PmuSelMtx = full(sparse(1:P.MaxNumPmu, P.PMU_POS, 1, P.MaxNumPmu, G.n));

% solve power flow for nominal case
[results, success] = runpf(G.mpc, P.mpopt);

% get nominal values
G.s = -(results.bus(:,PD) + 1i.*results.bus(:,QD))./results.baseMVA;
G.u = results.bus(:,VM).*exp(1i.*deg2rad(results.bus(:,VA)));

% currents
G.i = G.L * G.u;

% matrix for linearization
P.N = kron(eye(G.n),diag([1,-1]));
P.Ru = [diag(cos(angle(G.u))) , -diag(abs(G.u))*diag(sin(angle(G.u))) ; ...
      diag(sin(angle(G.u))) ,  diag(abs(G.u))*diag(cos(angle(G.u)))];
G.A_x_star = (c2ri(diag(conj(G.L*G.u))) + c2ri(diag(G.u))*P.N*c2ri(G.L))*P.Ru;
G.Apq = G.A_x_star(logical(kron(G.idx_pq,ones(2,1))), logical(kron(G.idx_pq,ones(2,1))));
G.Bpq = eye(size(G.Apq))/G.Apq;


%% generate priors, initial conditions in real-imag, true profiles and measurements
[Prior.s, Prior.u, Prior.i] = generate_priors(P, G);
% convert from complex to rectangular coordinates
IC.ds0 = zeros(2*G.num_pq, P.NumMC); 
for mc = 1:P.NumMC
    IC.ds0(:, mc) = c2ri(Prior.s(G.idx_pq, :, mc) - G.s(G.idx_pq), 'v');
end
IC.off0 = zeros(2*P.MaxNumPmu, 1);
IC.skew0 = zeros(2*P.MaxNumPmu, 1);

% true profiles
[True.s, True.u, True.i, True.off, True.skew] = generate_true_profiles(P, G);

% data and measurements
[Meas.delta, Meas.delay, Meas.u_nd, Meas.y_nd, Meas.u, Meas.y] = generate_measurements(P, G, True);


%% estimation

% Standard Estimation
[SE.s, SE.off, SE.skew, SE.u, SE.P, SE.TR.V, SE.ERR.V, SE.TR.OFF, ...
    SE.ERR.OFF, SE.TR.SKEW, SE.ERR.SKEW] = estimation_SASE(P, G, True, Meas, IC);
        
% Estimation assuming perfect knowledge of delay
[PFK.s, PFK.u, PFK.P, PFK.TR.V, PFK.ERR.V] = ...
     estimation_with_perfect_delay_knowledge(P, G, True, Meas, IC);
 
% Estimation assuming no knowledge of delay
[ABN.s, ABN.u, ABN.P, ABN.TR.V, ABN.ERR.V] = ...
   estimation_with_no_delay_knowledge(P, G, True, Meas, IC);

% Theoretical Posterior Values
if strcmp(exe_mode, mode_1)
    [Sigma, Std] = computeTheoreticalParametersPosterior(P, Meas.delta);
end                      


%% save
if saveData
    switch matpower_case
        case 'case15'
            save(strcat('data/CASE15_', date, save_suffix, '.mat'),'-v7.3');
        case 'case123'
            save(strcat('data/IEEE123_', date, save_suffix, '.mat'),'-v7.3');
    end
end

%% plot
% uncomment for inplace plotting or run separately
switch exe_mode
    case mode_1
        plot_IEEE_FuncNumPmuMeasurements;
    
    case mode_2
        plot_IEEE_FuncNumPmuDeployed;
    
    otherwise
        error('Selected exe_mode not supported');
end
