function [x, P] = KalmanFilterSASE(F, Q, H, R, y, x0, P0)
% Synchronization Aware State Estimator
% 
% Args:
%   F, Q, H, R: matrix for the Kalman Filter
%   y: measurements
%   x0, P0: initial condition and covariance
%
% Returns:
%   x, P: estimate of state and error covariance matrix

StateDim = length(x0);
[~,NumIter] = size(y);

if size(F,3) == 1
    F = repmat(F, [1, 1, NumIter]);
end

if size(R,3) == 1
    R = repmat(R, [1, 1, NumIter]);
end

if size(H,3) == 1
    H = repmat(H, [1, 1, NumIter]);
end

if size(Q,3) == 1
    Q = repmat(Q, [1, 1, NumIter]);
end

% support matrices
I = eye(StateDim);

% init
xt = sparse(x0);
Pt = sparse(P0);

% memories
x = zeros(NumIter, StateDim);
P = zeros(StateDim, StateDim, NumIter);

% for loop
for t = 1:NumIter
        
    % current matrices
    Ft = sparse(F(:,:,t));
    Qt = sparse(Q(:,:,t));
    Ht = sparse(H(:,:,t));
    Rt = sparse(R(:,:,t));
    yt = sparse(y(:,t));
            
    % prediction
    xt = sparse(Ft*xt);
    Pt = sparse(Ft*Pt*Ft' + Qt);
    
    % correction (no prediction given static model and no input)
    K = (Pt*Ht')/(Ht*Pt*Ht' + Rt);
    xt = xt + K*(yt - Ht*xt);
    Pt = (I - K*Ht)*Pt*(I - K*Ht)' + K*Rt*K';

    % save in complex coordinate
    x(t,:) = xt;
    P(:,:,t) = Pt;
    
end

end

