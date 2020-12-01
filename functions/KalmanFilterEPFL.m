function [x_abs, x_ang, P] = KalmanFilterEPFL(F, H, R, y, x0, P0)
% Run KF according to M. Pignati, et al., ISGT2015
%   Real-Time State Estimation of the EPFL-Campus Medium-Voltage Grid by Using PMUs. 
%   The Sixth Conference on Innovative Smart Grid Technologies (ISGT2015), 2015.
%
% Args:
%   F, H, R: state space and mesurement matrix for KF model
%   y: measurements
%   x0: initial condition
%   P0: initial state covariance
%
% Returns:
%   x_abs, x_ang: estimate in abs and phase
%   P: covariance matrix

[~, NumIter] = size(y);
NumNodes = size(F, 1) / 2;

if size(R, 3) == 1
    R = repmat(R, [1, 1, NumIter]);
end
if size(F, 3) == 1
    F = repmat(F, [1, 1, NumIter]);
end

% initialization
xt = x0;
Pt = P0;

% memories
x = zeros(NumIter, NumNodes);
P = zeros(2 * NumNodes, 2 * NumNodes, NumIter); 

% for loop over the measurements
In = eye(2 * NumNodes);
Q = zeros(2 * NumNodes);
xT = zeros(50, 2 * NumNodes);
for t = 1:NumIter
    
    % measurement and state matrix
    Ft = F(:, :, t);
    Rt = R(:, :, t);
    
    % update process covariance matrix
    if t > 50
        xOld = x(t-50 :t -1,:);
        xT(:, 1:2:2*NumNodes) = real(xOld);
        xT(:, 2:2:2*NumNodes) = imag(xOld);
        Q = diag(diag(cov(xT)));
    end
    
    % prediction
    xt = Ft*xt;
    Pt = Ft*Pt*Ft.' + Q;

    % correction (no prediction given static model and no input)
    K = (Pt*H')*pinv(H*Pt*H' + Rt);
    xt = xt + K*(y(:,t) - H*xt);
    Pt = (In - K*H)*Pt*(In - K*H)' + K*Rt*K';
    
    % save in complex coordinate
    x(t,:) = xt(1:2:2*NumNodes) + 1i*xt(2:2:2*NumNodes);
    P(:,:,t) = Pt;
    
end

x_abs = abs(x);
x_ang = unwrap(angle(x));

end

