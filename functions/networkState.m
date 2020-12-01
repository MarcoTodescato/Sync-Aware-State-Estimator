function [u, i, success] = networkState(G, s, sc)
% [u, i] = networkState(G,s,sc)
%  computes the exact solution of the network power flow analysis for
%  the grid G, where s is the vector of nominal complex power injected by
%  the nodes, and sc is the vector of complex power injected by the 
%  compensators.
%
%  Args:
%   G descibes the grid
%   s is the complex power injected by each node
%   sc is the complex power injected by compensators
%
%  Returns:  
%   u is the voltage at each node
%   i is the current injected by each node
%
%   s, sc, u, i have length G.n

% initialize at the approximate network state
u = ones(G.n, 1) * G.U0;

% threshold
MAX_REL_POWER_ERR = 1e-10;

epsil = 0.5;
i = conj((s+sc) ./ u);

% error
er = Inf;

% flag
success = 1;

MAX_ITER = 1e4;
iter = 0;

while (max(er) > MAX_REL_POWER_ERR) && (iter < MAX_ITER) 

    % update num of iter
    iter = iter + 1;
    
    % complex power at this voltage
    target_s = s.*(abs(u)/abs(G.U0)).^G.eta + sc;

    % current that gives target_s at this voltage
    i = (1-epsil)*i + epsil*conj(target_s ./ u);

    % voltage with the new current injection
    u = G.U0 + G.X * i;

    % error
    er =  target_s - u .* conj(i);
    er(G.PCC) = 0;
    er = er ./ target_s;

end

if iter >= MAX_ITER
    success = 0;
end


