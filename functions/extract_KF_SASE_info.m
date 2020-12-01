function [v, theta, v_pf, theta_pf, offset, skew, Pv, Ptheta, Poff, Pskew,...
          p, q, Pp, Pq] = extract_KF_SASE_info(LblPMU, Lbl0Inj, G, V, I, S, x, P)
% From the estimated state and covariance matrices output of the 
% KalmanFilterSASE in Power coordinate, the functions returns the estimated 
% quantities in voltage coordinates with corresponding covariances.
% 
% Args:
%   LblPMU: label of nodes where PMU meas are taken from
%   Lbl0Inj: label of zero-injection node
%   G: struct with grid
%   V: complex voltage measurement
%   I: complex current measurement
%   S: complex power measurements
%   x, P: estimated state and error covariance matrix
%
% Returns:
%   v: estimated votage abs
%   theta: estimated voltage phase
%   v_pf: voltage abs computed via exact power flow
%   theta_pf: voltage phase via exact power flow
%   offset: estimated offset sync error parameter (beta)
%   skew: estimated skew sync error parameter (alpha)
%   Pv: abs voltage covariance 
%   Ptheta: phase voltage covariance
%   Poff: offset covariance
%   Pskew: skew covariance
%   p: estimated real power around nominal point
%   q: estimated imag power around nominal point
%   Pp: real power covariance
%   Pq: imag power covariance

PCC = G.PCC;
L = G.L;

TotNumNodes = size(L,1);
LblNo0Inj = setdiff(1:TotNumNodes,Lbl0Inj);
Lno0inj = L(LblNo0Inj,LblNo0Inj) - L(LblNo0Inj,Lbl0Inj)*(L(Lbl0Inj,Lbl0Inj)\L(Lbl0Inj,LblNo0Inj));

LblNo0Inj2 = sort([2*LblNo0Inj-1 2*LblNo0Inj]);

NumNodes = size(Lno0inj,1);
NumIter = size(V,1);
NumPMU = length(LblPMU);
ds = x(:,1:2*NumNodes)';
offset = x(:,2*NumNodes+(1:NumPMU));
skew = x(:,2*NumNodes+NumPMU+(1:NumPMU));

% normalize measurements wrt PCC
Vpcc = mean(abs(V(:,PCC))); % pcc mean absolute value
V = V./Vpcc;
S = S./Vpcc;    % scale powers
L = L.*Vpcc;    % scale laplacia matrix
Lno0inj = Lno0inj.*Vpcc;

% compute average/predicited/nominal working point
s0 = [mean(S,1).' ;0];
[u0,i0] = networkState(G,s0,zeros(TotNumNodes,1));
% u0 = u0(LblNo0Inj);
% i0 = i0(LblNo0Inj);

% linearization matrix
Id = eye(NumNodes);
Ru = ri2p_jacobian(u0(LblNo0Inj));
N = kron(Id,[1 0 ; 0 -1]);
Hv2s = (c2ri(diag(conj(i0(LblNo0Inj)))) + c2ri(diag(u0(LblNo0Inj)))*N*c2ri(Lno0inj))*Ru;
Hs2v = inv(Hv2s);

% estimated power
s_ri = repmat(c2ri(s0(LblNo0Inj),'v'),[1,NumIter]) + ds;
s = zeros(NumIter,TotNumNodes);
s(:,LblNo0Inj) = (s_ri(1:2:end,:) + 1i*s_ri(2:2:end,:)).';

% estimated voltage
dv = Hs2v*ds;
u = zeros(NumIter,TotNumNodes);
u(:,LblNo0Inj) = p2c(repmat(c2p(u0(LblNo0Inj)),[1,NumIter]) + dv).';
u(:,Lbl0Inj) = -u(:,LblNo0Inj)*((L(Lbl0Inj,Lbl0Inj)\L(Lbl0Inj,LblNo0Inj)).');

% covariances
Ps = nan(2*TotNumNodes,2*TotNumNodes,NumIter);
Ps(LblNo0Inj2,LblNo0Inj2,:) = P(1:2*NumNodes,1:2*NumNodes,:);

Pu = nan(2*TotNumNodes,2*TotNumNodes,NumIter);
for t = 1:NumIter
    Pu(LblNo0Inj2,LblNo0Inj2,t) = Hs2v*Ps(LblNo0Inj2,LblNo0Inj2,t)*Hs2v';
end

Poff = P(2*NumNodes+(1:NumPMU),2*NumNodes+(1:NumPMU),:);
Pskew = P(2*NumNodes+NumPMU+(1:NumPMU),...
          2*NumNodes+NumPMU+(1:NumPMU),:);
      
% de-normalize all the quantities
s = s.*Vpcc;
u = u.*Vpcc;
u_pf = zeros(NumIter,TotNumNodes);
for t = 1:NumIter
    [u_tmp,~] = networkState(G,s0,zeros(TotNumNodes,1));
    u_pf(t,:) = (u_tmp).';  
end

Ps = Ps*Vpcc^2;
Pu = Pu*Vpcc^2;

% extract desired quantities
p = real(s);
q = imag(s);
v = abs(u);
theta = unwrap(angle(u));
v_pf = abs(u_pf);
theta_pf = unwrap(angle(u_pf));
Pp = Ps(1:2:end,1:2:end,:);
Pq = Ps(2:2:end,2:2:end,:);
Pv = Pu(1:2:end,1:2:end,:);
Ptheta = Pu(2:2:end,2:2:end,:);
