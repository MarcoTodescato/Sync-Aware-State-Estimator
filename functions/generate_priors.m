function [s, u, i] = generate_priors(Param, Grid, varargin)
% Given a set of parameters and a testcase, generate prior info
%
% Args:
%   Param: struct of parameters
%   Grid: struct of testcase
%
% Returns:
%   s: complex powers
%   u: complex voltages
%   i: complex currents

% matpower constants
PD = 3;
QD = 4;
VM = 8;
VA = 9;

if nargin > 2
    PD = varargin{1};
    QD = varargin{2};
    VM = varargin{3};
    VA = varargin{4};
end

% init containers
s = zeros(Grid.n,1,Param.NumMC);
u = zeros(Grid.n,1,Param.NumMC);
i = zeros(Grid.n,1,Param.NumMC);

% matpower case
mpc = Grid.mpc;

% loop over the number of MonteCarlo runs
for mc = 1:Param.NumMC
    
    % generate powers for each monte carlo run from the nominal value of
    % the given test case
    s(:,:,mc) = Grid.s + (Param.real_prior_load_variation .* abs(real(Grid.s)) .* randn(Grid.n,1) + ...
                           1i*Param.imag_prior_load_variation .* abs(imag(Grid.s)) .* randn(Grid.n,1));
    
    % solve the load flow...
    switch Param.solver_type
        
        %... either with a linear approximation
        case 'linear'
            u(:,:,mc) = Grid.v0 + Grid.X*conj(s(:,:,mc))/Grid.v0;

        %... or exactly
        otherwise
            mpc.bus(:,PD) = -real(s(:,:,mc)).*mpc.baseMVA;
            mpc.bus(:,QD) = -imag(s(:,:,mc)).*mpc.baseMVA;
            [results, ~] = runpf(mpc,Param.mpopt);
            u(:,:,mc) = results.bus(:,VM).*exp(1i.*deg2rad(results.bus(:,VA)));
            
    end
    
    % get currents
    i(:,:,mc) = Grid.L*u(:,:,mc);
    
end

end

