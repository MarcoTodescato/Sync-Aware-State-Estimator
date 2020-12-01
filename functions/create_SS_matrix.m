function [A, B, C] = create_SS_matrix(Param, Grid, num_pmu, delta)
% Create SS matrices needed for SASE estimator
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

%A = zeros(2*(Grid.num_pq+2*num_pmu), 2*(Grid.num_pq+2*num_pmu), NumIter);
%B = zeros(2*(Grid.num_pq+2*num_pmu), 2*Grid.num_pq, NumIter);
C = zeros(2*num_pmu, 2*(Grid.num_pq+2*num_pmu), NumIter);
idx_delta = kron(ones(1,Param.NumWindow),1:length(delta));

switch Param.state_correlation
    case 'uncorrelated'
        At = c2ri(blkdiag(zeros(Grid.num_pq),eye(2*num_pmu)));
        Bt = c2ri([eye(Grid.num_pq) ; zeros(2*num_pmu,Grid.num_pq)]);
        
    case 'correlated'
        At = c2ri(blkdiag(eye(Grid.num_pq),eye(2*num_pmu)));
        Bt = c2ri([zeros(Grid.num_pq) ; zeros(2*num_pmu,Grid.num_pq)]);
end
A = repmat(At,1,1,NumIter);
B = repmat(Bt,1,1,NumIter);

idx = logical(kron(pad(Grid.num_pq,Param.PMU_POS(1:num_pmu) - 1,1),ones(2,1)));

%for t = 1:NumIter
    % state space model matrices
%     switch Param.state_correlation
%         case 'uncorrelated'
%             A(:,:,t) = c2ri(blkdiag(zeros(Grid.num_pq),eye(2*num_pmu)));
%             B(:,:,t) = c2ri([eye(Grid.num_pq) ; zeros(2*num_pmu,Grid.num_pq)]);
%             
%         case 'correlated'
%             A(:,:,t) = c2ri(blkdiag(eye(Grid.num_pq),eye(2*num_pmu)));
%             B(:,:,t) = c2ri([zeros(Grid.num_pq) ; zeros(2*num_pmu,Grid.num_pq)]);
%             
%     end

switch Param.coordinate_type
    case 'rectangular'
        for t = 1:NumIter
            
            C(:,:,t) = c2ri([X_pmu, 1i*eye(num_pmu) , ...
                1i*eye(num_pmu)*delta(idx_delta(t))]);
        end
    case 'polar'
        
        for t = 1:NumIter
            %idx = logical(kron(pad(Grid.num_pq,Param.PMU_POS(1:num_pmu) - 1,1),ones(2,1)));
            C(:,:,t) = [Grid.Bpq(idx,:), c2ri([1i*eye(num_pmu) , 1i*eye(num_pmu)*delta(t)])];
        end
        
end
%end

end

