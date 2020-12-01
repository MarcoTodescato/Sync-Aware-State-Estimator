function [S, U, EST_COV, TR_V, ERR_V] = ...
    estimation_with_perfect_delay_knowledge(Param, Grid, TrueProfiles, Meas, IC)
% Standard KF on a synthetic loadcase assuming measurements expressed in
% polar coordianates and perferct delay knowledge (which is compensanted)
% given parameters, loadcase, true profiles, measurements and initial conditions
%
% Args:
%   Param: struct of parameters
%   Grid: struct with loadcase
%   TrueProfiles: previously generated true electric profiles
%   Meas: previously generated measurements
%   IC: previously generated initial conditions
%
% Returns:
%   S: estimated power profiles
%   U: estimated voltage profiles
%   EST_COV: complete estimated state error covariance matrix
%   TR_V: trace of error covariance matrix for estimated voltage profiles
%   ERR_V: MSE between true and estimated voltage

% memories
S           = zeros(Grid.num_pq,length(Param.NUM_MEAS),Param.NumMC);
U           = zeros(Grid.num_pq,length(Param.NUM_MEAS),Param.NumMC);
EST_COV     = zeros(2*Grid.num_pq,2*Grid.num_pq,length(Param.NUM_MEAS));
TR_V        = zeros(length(Param.NUM_PMU),length(Param.NUM_MEAS),Param.NumMC);
ERR_V       = zeros(length(Param.NUM_PMU),length(Param.NUM_MEAS),Param.NumMC);

% loop over the number of MC
for mc_it = 1:Param.NumMC
    
    % loop over the number of PMU
    for pmu_it = 1:length(Param.NUM_PMU);
        
        num_pmu = Param.NUM_PMU(pmu_it);
        
        for meas_it = 1:length(Param.NUM_MEAS)
            
            % current number of measurement taken
            m = Param.NUM_MEAS(meas_it);
            NumIter = m*Param.NumWindow;
            
            % vector to extract desired measurements
            pos_delta = round(linspace(1,Param.MaxNumMeas,m))';
            pos_meas = kron(ones(Param.NumWindow,1),pos_delta) + ...
                kron((0:Param.NumWindow-1)',Param.MaxNumMeas*ones(m,1));
            
            % state space model matrices
            [A,B,C] = create_SS_matrix_reduced(Param,Grid,num_pmu,Meas.delta(pos_delta));
            
            % noise covariance matrix
            Cov_s = (diag(kron(ones(Grid.num_pq,1),...
                [Param.real_percent_load_variation;Param.imag_percent_load_variation]))*...
                diag(abs(c2ri(Grid.s(Grid.idx_pq),'v'))))^2;
            Cov_s0 = (diag(kron(ones(Grid.num_pq,1),...
                [Param.real_prior_load_variation;Param.imag_prior_load_variation]))*...
                diag(abs(c2ri(Grid.s(Grid.idx_pq),'v'))))^2;
            R = c2ri(((Param.sigma_pmu*Grid.v0)^2*eye(num_pmu) + ...
                Param.sigma_pmu^2*eye(num_pmu))/2);
            
            % build input output statistic
            Q = zeros(2*Grid.num_pq,2*Grid.num_pq,NumIter);
            Rt = zeros(2*num_pmu,2*num_pmu,NumIter);
            for t = 1:NumIter
                Q(:,:,t) = Cov_s;
                Rt(:,:,t) = R;
            end
            u = IC.ds0(:,mc_it)*ones(1,NumIter);
            y = Meas.y_nd(1:2*num_pmu,pos_meas,mc_it);
            
            %  priors for real and imaginary parts
            x0 = IC.ds0(:,mc_it);
            P0 = Cov_s0;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % kalman filtering
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x = zeros(Grid.num_pq,NumIter);
            P = zeros(2*Grid.num_pq,2*Grid.num_pq,NumIter);
            
            xt = x0;
            Pt = P0;
            In = eye(length(xt));
            for t = 1:NumIter
                
                % prediction
                xt = A(:,:,t)*xt + B(:,:,t)*u(:,t);
                Pt = A(:,:,t)*Pt*A(:,:,t)' + Q(:,:,t);
                
                % correction
                K = (Pt*C(:,:,t)')/(C(:,:,t)*Pt*C(:,:,t)' + Rt(:,:,t));
                xt = xt + K*(y(:,t) - C(:,:,t)*xt);
                Pt = (In - K*C(:,:,t))*Pt*(In - K*C(:,:,t))' + K*Rt(:,:,t)*K';
                
                % save in complex coordinate
                x(:,t) = kron(eye((Grid.num_pq)),[1 1i])*xt;
                P(:,:,t) = Pt;
            end
            ds_hat = x(:,end);
            u_hat = p2c(c2p(Grid.u(Grid.idx_pq)) + Grid.Bpq*c2ri(ds_hat,'v'));
            
            % errors
            TR_V(pmu_it,meas_it,mc_it)  = sqrt(trace(Grid.Bpq*P(:,:,end)*Grid.Bpq')./Grid.num_pq);
            
            ERR_V(pmu_it,meas_it,mc_it) = mean(abs(TrueProfiles.u(Grid.idx_pq,end,mc_it) - u_hat).^2);
                        
            % save estimation
            if num_pmu == Param.MaxNumPmu
                S(:,meas_it,mc_it) = ds_hat + TrueProfiles.s(Grid.idx_pq,end,mc_it);
                U(:,meas_it,mc_it) = u_hat;
                EST_COV(:,:,meas_it) = P(:,:,end);
            end
            
        end
        
    end
end

%TR_V  = mean(TR_V,2);
%ERR_V = sqrt(mean(ERR_V,2));
TR_V  = mean(TR_V,3);
ERR_V = sqrt(mean(ERR_V,3));


