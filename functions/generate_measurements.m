function [delta, delay, u_meas_nodelay, y_nd, u_meas, y] = ...
    generate_measurements(Param, Grid, TrueProfiles)
% Given a set of parameters and a testcase, and previously generated true
% electric profiles, generate mesurement data for estimation
%
% Args:
%   Param: struct of parameters
%   Grid: struct of testcase
%   TrueProfiles: struct with underlying true electric profiles
%   measurements are sampled from
%
% Returns:
%   delta: pmu sampling instants
%   delay: sync error
%   u_meas_nodelay: voltage measurements without delay error
%   y_nd: measurements w/o delay either in rectangular or polar coordinates
%       depending on selected parameter
%   u_meas: voltage measurements with delay
%   y: measurements with delay either in rectangular or polar coordinates
%       depending on selected parameter


% parameters for the PMU data
delta = Param.tau / Param.MaxNumMeas * (1:1:Param.MaxNumMeas)';

% generate delay
delay = zeros(Param.MaxNumPmu, Param.MaxNumIter, Param.NumMC);
for mc = 1:Param.NumMC
    for w = 1:Param.NumWindow
        delay(:,(w-1)*Param.MaxNumMeas+(1:Param.MaxNumMeas),mc) = ...
             TrueProfiles.off(:,w,mc)*ones(1,Param.MaxNumMeas) + ...
             TrueProfiles.skew(:,w,mc)*delta';
    end
end

% generate measurements
switch Param.meas_model
    case 'true'
            abs_noise = randn(Param.MaxNumPmu,Param.MaxNumIter,Param.NumMC);
            angle_noise = randn(Param.MaxNumPmu,Param.MaxNumIter,Param.NumMC);
            u_meas_nodelay = (abs(TrueProfiles.u(Param.PMU_POS,:,:)) + ...
                                       Param.sigma_pmu*Grid.v0*abs_noise).*...
                                    exp(1i.*(angle(TrueProfiles.u(Param.PMU_POS,:,:)) + ...
                                    Param.sigma_pmu*angle_noise));
            u_meas = (abs(TrueProfiles.u(Param.PMU_POS,:,:)) + ...
                                Param.sigma_pmu*Grid.v0*abs_noise).*...
                            exp(1i.*(angle(TrueProfiles.u(Param.PMU_POS,:,:)) + ...
                            Param.sigma_pmu*angle_noise + delay));    
        
    otherwise
            u_meas_nodelay = TrueProfiles.u(Param.PMU_POS,:,:) + ...
                     Param.sigma_pmu*(randn(Param.MaxNumPmu,Param.MaxNumIter,Param.NumMC)  + ...
                                    1i*randn(Param.MaxNumPmu,Param.MaxNumIter,Param.NumMC));
            u_meas =  u_meas_nodelay + 1i*delay;
end

% switch to phase-angle (polar) or rectangular
switch Param.coordinate_type
    
    case 'rectangular'
        y_nd = c2ri(u_meas_nodelay - ...
                Grid.v0*ones(Param.MaxNumPmu,Param.MaxNumIter,Param.NumMC),'evol');
        y = c2ri(u_meas - Grid.v0*ones(Param.MaxNumPmu,Param.MaxNumIter,Param.NumMC),'evol');
        
    case 'polar'
        y_nd = c2p(u_meas_nodelay) - superkron(c2p(Grid.u(Param.PMU_POS)),ones(1,Param.MaxNumMeas,Param.NumMC));
        y = c2p(u_meas) - superkron(c2p(Grid.u(Param.PMU_POS)),ones(1,Param.MaxNumMeas,Param.NumMC));
end

end

