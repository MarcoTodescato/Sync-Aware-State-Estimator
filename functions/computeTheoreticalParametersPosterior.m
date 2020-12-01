function [Sigma, Std] = computeTheoreticalParametersPosterior(Parameters, delta)
% Given noise parameters, compute theoretical covariance and std values

% init the memory
Sigma = zeros(3,3,length(Parameters.NUM_MEAS));
Std = zeros(3,length(Parameters.NUM_MEAS));

% create prior variance
B = [1 -1 0; 1 1 0; 0 0 1];
P0 = sparse(B*diag([Parameters.sigma_off(1),...
                    Parameters.real_prior_load_variation(1),...
                    Parameters.sigma_skew(1)].^2)*B');

% create prior information
I0 = inv(P0);

% noise variance matrix
R = Parameters.sigma_pmu^2;

for meas_it = 1:length(Parameters.NUM_MEAS)
    
    M = Parameters.NUM_MEAS(meas_it);
    
    % number of measurements
    pos_delta = round(linspace(1,Parameters.MaxNumMeas,M))';
    
    % measurements matrix
    C = sparse([zeros(M,1) ones(M,1) delta(pos_delta)]);

    % measurement info
    Im = C'*(R\C);
    
    % posterior info
    I = I0 + Im;
    
    % posterior variance
    Sigma(:,:,meas_it) = B\(eye(3)/I)/B';
    Std(:,meas_it) = diag(sqrtm(Sigma(:,:,meas_it)));

end

