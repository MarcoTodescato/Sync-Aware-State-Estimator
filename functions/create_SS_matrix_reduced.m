function [A, B, C] = create_SS_matrix_reduced(Param, Grid, num_pmu, delta)
% Create SS matrices needed for standard KF (no sync error parameter
% estimation)
%
% Args:
%   Param: struct of param
%   Grid: struct with loadcase
%   num_pmu: number of deployed pmu
%   delta: pmu sampling instant
%
% Returns:
%   A, B, C: state, input and output matrices

pos_pmu = Param.PMU_POS(1:num_pmu);
X_pmu = Grid.X(pos_pmu,Grid.idx_pq);
NumIter = length(delta)*Param.NumWindow;

A = zeros(2*Grid.num_pq, 2*Grid.num_pq, NumIter);
B = zeros(2*Grid.num_pq, 2*Grid.num_pq, NumIter);
C = zeros(2*num_pmu,  2*Grid.num_pq, NumIter);

for t = 1:NumIter
    % state space model matrices
    
    switch Param.state_correlation
        case 'uncorrelated'
            A(:,:,t) = c2ri(zeros(Grid.num_pq));
            B(:,:,t) = c2ri(eye(Grid.num_pq));
            
        case 'correlated'
            A(:,:,t) = c2ri(eye(Grid.num_pq));
            B(:,:,t) = c2ri(zeros(Grid.num_pq));
    end
    
    switch Param.coordinate_type
        case 'rectangular'
            C(:,:,t) = c2ri(X_pmu);
        
        case 'polar'
            idx = logical(kron(pad(Grid.num_pq,Param.PMU_POS(1:num_pmu) - 1,1),ones(2,1)));
            C(:,:,t) = Grid.Bpq(idx,:);
    end
end

end

