function [S, OFF, SKEW, U, EST_COV, TR_V, ERR_V, TR_OFF, ERR_OFF, TR_SKEW, ERR_SKEW] = ....
            estimation_SASE(Param, Grid, TrueProfiles, Meas, IC)
% Proposed SASE KF on a synthetic loadcase assuming measurements expressed in polar coordianates 
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
%   OFF: estimated offset parameter (beta)
%   U: estimated voltage profiles
%   EST_COV: complete estimated state error covariance matrix
%   TR_V: trace of error covariance matrix for estimated voltage profiles
%   ERR_V: MSE between true and estimated voltage
%   TR_OFF: trace of error covariance matrix for estimated offset profiles
%   ERR_OFF: MSE between true and estimated offset
%   TR_SKEW: trace of error covariance matrix for estimated skew profiles
%   ERR_SKEW: MSE between true and estimated skew

% memories
S           = zeros(Grid.num_pq,length(Param.NUM_MEAS),Param.NumMC);
OFF         = zeros(Param.MaxNumPmu,length(Param.NUM_MEAS),Param.NumMC);
SKEW        = zeros(Param.MaxNumPmu,length(Param.NUM_MEAS),Param.NumMC);
U           = zeros(Grid.num_pq,length(Param.NUM_MEAS),Param.NumMC);
EST_COV     = zeros(2*Grid.num_pq+4*Param.MaxNumPmu,2*Grid.num_pq+4*Param.MaxNumPmu,...
                    length(Param.NUM_MEAS));
TR_OFF      = zeros(length(Param.NUM_PMU),length(Param.NUM_MEAS),Param.NumMC);
TR_SKEW     = zeros(length(Param.NUM_PMU),length(Param.NUM_MEAS),Param.NumMC);
TR_V        = zeros(length(Param.NUM_PMU),length(Param.NUM_MEAS),Param.NumMC);
ERR_OFF     = zeros(length(Param.NUM_PMU),length(Param.NUM_MEAS),Param.NumMC);
ERR_SKEW    = zeros(length(Param.NUM_PMU),length(Param.NUM_MEAS),Param.NumMC);
ERR_V       = zeros(length(Param.NUM_PMU),length(Param.NUM_MEAS),Param.NumMC);


% loop over the number of MC
for mc_it = 1:Param.NumMC
    
    for pmu_it = 1:length(Param.NUM_PMU);
    
        num_pmu = Param.NUM_PMU(pmu_it);
        
        % loop over the num of meas per PMU
        for meas_it = 1:length(Param.NUM_MEAS)
            
            % current number of measurement taken
            m = Param.NUM_MEAS(meas_it);
            NumIter = m*Param.NumWindow; 
            
            % vector to extract desired measurements
            pos_delta = round(linspace(1,Param.MaxNumMeas,m))';
            pos_meas = kron(ones(Param.NumWindow,1),pos_delta) + ...
                        kron((0:Param.NumWindow-1)',Param.MaxNumMeas*ones(m,1));
            
            % state space model matrices
            [A,B,C] = create_SS_matrix(Param,Grid,num_pmu,Meas.delta(pos_delta));
    
            % noise covariance matrix
            Cov_off  = kron(diag(Param.sigma_off^2.*ones(num_pmu,1)),diag([1;0]));
            Cov_skew = kron(diag(Param.sigma_skew^2.*ones(num_pmu,1)),diag([1;0]));
            Cov_s = (diag(kron(ones(Grid.num_pq,1),...
                          [Param.real_percent_load_variation;Param.imag_percent_load_variation]))*...
                diag(abs(c2ri(Grid.s(Grid.idx_pq),'v'))))^2;
            Cov_s0 = (diag(kron(ones(Grid.num_pq,1),...
                            [Param.real_prior_load_variation;Param.imag_prior_load_variation]))*...
                diag(abs(c2ri(Grid.s(Grid.idx_pq),'v'))))^2;
            R = c2ri(((Param.sigma_pmu*Grid.v0)^2*eye(num_pmu) + ...
                Param.sigma_pmu^2*eye(num_pmu))/2);
            
            % build input output statistic           
            Q = zeros(2*(Grid.num_pq)+4*num_pmu,2*(Grid.num_pq)+4*num_pmu,NumIter);
            Rt = zeros(2*num_pmu,2*num_pmu,NumIter);
            for t = 1:NumIter
                Q(:,:,t) = blkdiag(Cov_s,c2ri(zeros(2*num_pmu)));
                Rt(:,:,t) = R;
            end
            u = IC.ds0(:,mc_it)*ones(1,NumIter);
            y = Meas.y(1:2*num_pmu,pos_meas,mc_it);
           
            %  priors for real and imaginary parts
            x0 = [IC.ds0(:,mc_it); IC.off0(1:2*num_pmu) ; IC.skew0(1:2*num_pmu)];
            P0 = blkdiag(Cov_s0,Cov_off,Cov_skew);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % kalman filtering
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x = zeros((Grid.num_pq)+2*num_pmu,NumIter);
            P = zeros(2*(Grid.num_pq)+4*num_pmu,2*(Grid.num_pq)+4*num_pmu,NumIter);
            
            xt = x0;
            Pt = P0;
            Pt_tmp = P0;
            In = eye(length(xt));
            for t = 1:NumIter

                Pt_tmp = A(:,:,t)*Pt_tmp*A(:,:,t)';
                K_tmp = (Pt_tmp*C(:,:,t)')/(C(:,:,t)*Pt_tmp*C(:,:,t)' + Rt(:,:,t));
                Pt_tmp = (In - K_tmp*C(:,:,t))*Pt_tmp*(In - K_tmp*C(:,:,t))' + K_tmp*Rt(:,:,t)*K_tmp';
                
                % prediction
                xt = A(:,:,t)*xt + B(:,:,t)*u(:,t);
                Pt = A(:,:,t)*Pt*A(:,:,t)' + Q(:,:,t);
                
                % correction
                K = (Pt*C(:,:,t)')/(C(:,:,t)*Pt*C(:,:,t)' + Rt(:,:,t));
                xt = xt + K*(y(:,t) - C(:,:,t)*xt);
                Pt = (In - K*C(:,:,t))*Pt*(In - K*C(:,:,t))' + K*Rt(:,:,t)*K';
                    
                % save in complex coordinate
                x(:,t) = kron(eye((Grid.num_pq)+2*num_pmu),[1 1i])*xt;
                P(:,:,t) = Pt;
            end
            ds_hat = x(1:Grid.num_pq,end);
            off_hat = x(Grid.num_pq+1:Grid.num_pq+num_pmu,end);
            skew_hat = x(Grid.num_pq+num_pmu+1:end,end);
            u_hat = p2c(c2p(Grid.u(Grid.idx_pq)) + Grid.Bpq*c2ri(ds_hat,'v'));            

            % errors
            TR_OFF(pmu_it,meas_it,mc_it)  = sqrt(trace(P(2*Grid.num_pq+(1:2*num_pmu),...
                                                        2*Grid.num_pq+(1:2*num_pmu),end))./num_pmu);
            TR_SKEW(pmu_it,meas_it,mc_it) = sqrt(trace(P(2*Grid.num_pq+2*num_pmu+1:end,...
                                                        2*Grid.num_pq+2*num_pmu+1:end,end))./num_pmu);
            TR_V(pmu_it,meas_it,mc_it)    = sqrt(trace(Grid.Bpq*...
                                                          P(1:2*Grid.num_pq,1:2*Grid.num_pq,end)*...
                                                          Grid.Bpq')./Grid.num_pq);
                        
            ERR_OFF(pmu_it,meas_it,mc_it)  = mean(abs(TrueProfiles.off(1:num_pmu,end,mc_it) - off_hat).^2);
            ERR_SKEW(pmu_it,meas_it,mc_it) = mean(abs(TrueProfiles.skew(1:num_pmu,end,mc_it) - skew_hat).^2);            
            ERR_V(pmu_it,meas_it,mc_it)    = mean(abs(TrueProfiles.u(Grid.idx_pq,end,mc_it) - u_hat).^2);
                        
            % save estimation
            if num_pmu == Param.MaxNumPmu
                S(:,meas_it,mc_it) = ds_hat + TrueProfiles.s(Grid.idx_pq,end,mc_it);
                OFF(:,meas_it,mc_it) = off_hat;
                SKEW(:,meas_it,mc_it) = skew_hat;
                U(:,meas_it,mc_it) = u_hat;
                EST_COV(:,:,meas_it) = P(:,:,end); 
            end
        end
    end
end

TR_OFF      = mean(TR_OFF,3);
TR_SKEW     = mean(TR_SKEW,3);
TR_V        = mean(TR_V,3);
ERR_OFF     = sqrt(mean(ERR_OFF,3));
ERR_SKEW    = sqrt(mean(ERR_SKEW,3));
ERR_V       = sqrt(mean(ERR_V,3));
