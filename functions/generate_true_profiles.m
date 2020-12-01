function [s, u, i, offset, skew] = generate_true_profiles(Param, Grid, varargin)
% Given a set of parameters and a testcase, generate electric profiles and
% sync error parameters profiles used as underlying true profiles during estimation 
%
% Args:
%   Param: struct of parameters
%   Grid: struct of testcase
%
% Returns:
%   s: complex powers
%   u: complex voltages
%   i: complex currents
%   offset: sync error offset parameter (beta)
%   skew: sync error skew parameter (alpha)


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
s = zeros(Grid.n,Param.MaxNumIter,Param.NumMC);
u = zeros(Grid.n,Param.MaxNumIter,Param.NumMC);
i = zeros(Grid.n,Param.MaxNumIter,Param.NumMC);

% load case
mpc = Grid.mpc;

% pick state space model
switch Param.state_correlation
    case 'correlated'
        A = eye(Grid.n);
        B = zeros(Grid.n);
        
    otherwise
        A = zeros(Grid.n);
        B = eye(Grid.n);
end

% generate load profiles
for mc = 1:Param.NumMC
    % first iteration
    s(:,1,mc) = Grid.s + (Param.real_percent_load_variation.*abs(real(Grid.s)).*randn(Grid.n,1) + ...
                        1i*Param.imag_percent_load_variation.*abs(imag(Grid.s)).*randn(Grid.n,1));
    
    % depending on the selected type of SS model (matrices A and B) the 
    % following values will be (un)correlated
    for t = 2:Param.MaxNumIter
        s(:,t,mc) = A*s(:,t-1,mc)  + B*Grid.s + ...
                     (Param.real_percent_load_variation.*abs(real(Grid.s)).*randn(Grid.n,1) + ...
                      1i*Param.imag_percent_load_variation.*abs(imag(Grid.s)).*randn(Grid.n,1));
    end
end

% generate corresponding voltage and currents
switch Param.solver_type
    
    case 'linear'

        for mc = 1:Param.NumMC
            for t = 1:Param.MaxNumIter

                u(:,t,mc) = Grid.v0 + Grid.X*conj(s(:,t,mc))/Grid.v0;
                i(:,:,mc) = Grid.L*u(:,:,mc);
            end
        end

    otherwise
        
        for mc = 1:Param.NumMC
            for t = 1:Param.MaxNumIter
                mpc.bus(:,PD) = -real(s(:,t,mc)).*mpc.baseMVA;
                mpc.bus(:,QD) = -imag(s(:,t,mc)).*mpc.baseMVA;
                [results, ~] = runpf(mpc,Param.mpopt);
                u(:,t,mc) = results.bus(:,VM).*exp(1i.*deg2rad(results.bus(:,VA)));
                i(:,t,mc) = Grid.L*u(:,t,mc);
            end
        end
        
end
    

% generate offset and skew
%offset = Param.sigma_off*randn(Param.MaxNumPmu,Param.NumWindow,Param.NumMC);
%skew = Param.sigma_skew*randn(Param.MaxNumPmu,Param.NumWindow,Param.NumMC);
offset = zeros(Param.MaxNumPmu,Param.NumWindow,Param.NumMC);
skew = zeros(Param.MaxNumPmu,Param.NumWindow,Param.NumMC);
for w = 1:Param.NumWindow
    for mc = 1:Param.NumMC
        offset(:, w, mc) = Param.sigma_off .* randn(Param.MaxNumPmu, 1);
        skew(:, w, mc) = Param.sigma_skew .* randn(Param.MaxNumPmu, 1);
    end
end

end

