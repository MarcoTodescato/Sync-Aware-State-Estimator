function R = computeCovMatrixEPFLMeasurementsSASE(LabelPMU, PCC, baseV, baseA, V, I)
% Compute the covariance matrix for the EPFL measurements
% conveniently expressed in a form suitable to SASE
% 
% Args:
%   LabelPMU: label of nodes where PMU measurements are taken from
%   PCC: label of PCC node
%   baseV: base voltage value
%   baseA: base current value
%   V: complex voltage measurements
%   I: complex current measurements
%
% Returns:
%   R: covariance matrix


NumPMU = length(LabelPMU);
NumMeas = 4*NumPMU;
NumIter = size(V,1);

R = zeros(NumMeas,NumMeas,NumIter);


for t = 1:NumIter
    
    % complex voltage and current measurements ...
    ut = V(t,LabelPMU).';
    it = I(t,LabelPMU).';
    
    % ...into real and imaginary parts
    vre = real(ut);
    vim = imag(ut);
    
    ire = real(it);
    iim = imag(it);
    
    % build measurement covariance matrix
    r_v = zeros(2*NumPMU,1);
    r_s = zeros(2*NumPMU,1);
    
    [std_vm, std_va, std_im, std_ia, std_v_re, std_v_im, std_i_re, ...
        std_i_im] = computeStdMeasEPFL(ut,it,LabelPMU,PCC,baseV,baseA);
    
    r_v(1:2:end) = std_vm.^2;
    r_v(2:2:end) = std_va.^2;
    Rv = diag(r_v);
    
    r_s(1:2:end)  = vre.^2.*std_i_re.^2 + vim.^2.*std_i_im.^2 + ...
        ire.^2.*std_v_re.^2 + iim.^2.*std_v_im.^2;
    r_s(2:2:end)  = vre.^2.*std_i_im.^2 + vim.^2.*std_i_re.^2 + ...
        ire.^2.*std_v_im.^2 + iim.^2.*std_v_re.^2;
    
    %cross correlation among powers
    r_pq = vre.*vim.*(std_i_re.^2 - std_i_im.^2) + ...
        ire.*iim.*(std_v_re.^2 - std_v_im.^2);
    Rs = diag(r_s);
    for n = 1:NumPMU
        Rs(2*n-1,2*n) = r_pq(n);
        Rs(2*n,2*n-1) = r_pq(n);
    end
    
    R(:,:,t) = blkdiag(Rv, Rs);
end

end

