function [std_vm, std_va, std_im, std_ia, ...
          std_v_re, std_v_im, std_i_re, std_i_im] = computeStdMeasEPFL(v, i, bus_label, PCC, baseV, baseA)
% Computes the standard deviations to be used to compute the KF measurements covariance matrix
% according to
%   M. Pignati, et al. Real-Time State Estimation of the EPFL-Campus Medium-Voltage Grid by Using PMUs. 
%   The Sixth Conference on Innovative Smart Grid Technologies (ISGT2015), 2015.
%
% Args:
%   v: complex voltage
%   i: complex current
%   bus_label: label of buses for which std must be computed
%   PCC: index of PCC node
%   baseV: value for base voltage
%   baseA: value for base current
%
% Returns:
%   std for different electric quantities 


% standard deviation for volage
NumMeasVoltage = length(v);

vm = abs(v);
va = unwrap(angle(v));

std_vm = zeros(NumMeasVoltage,1);
std_va = zeros(NumMeasVoltage,1);

std_v_re = zeros(NumMeasVoltage,1);
std_v_im = zeros(NumMeasVoltage,1);

for m = 1:NumMeasVoltage
    
    rateV = 20000/(baseV); % this, in non normalized voltage, is 20KV
    
    % Linear Piece wise definition of voltage magnitude and phase meas std
    x_v = vm(m)/rateV*100;
    if x_v < 0.2
        f_vm = 10;
        f_va = 100;
        
    else if  x_v >= 0.2 && x_v < 5
            f_vm = -5/4.8*x_v + 10.2083;
            f_va = -70/4.8*x_v + 102.9167;
            
        else if x_v >=5 && x_v < 80
                f_vm = -4.9/75*x_v + 5.3267;
                f_va = -28.5/75*x_v + 31.9;
                
            else if x_v >= 80
                    f_vm = 0.1;
                    f_va = 1.5;
                end
            end
        end
    end
    
    % std for abs and phase
    std_vm(m) = (f_vm / 100) * vm(m) / 3;
    std_va(m) = f_va / 3 * 1e-3;
            
    % std for real and imag
    std_v_re(m) = sqrt(...
                   vm(m)^2*exp(-2*std_va(m).^2)*...
                        (cos(va(m))^2*(cosh(2*std_va(m)^2) - cosh(std_va(m)^2)) + ...
                         sin(va(m))^2*(sinh(2*std_va(m)^2) - sinh(std_va(m)^2))) + ...
                   std_vm(m)^2*exp(-2*std_va(m).^2)*...
                        (cos(va(m))^2*(2*cosh(2*std_va(m)^2) - cosh(std_va(m)^2)) + ...
                         sin(va(m))^2*(2*sinh(2*std_va(m)^2) - sinh(std_va(m)^2))));
    
   std_v_im(m) = sqrt(...
                  vm(m)^2*exp(-2*std_va(m).^2)*...
                        (sin(va(m))^2*(cosh(2*std_va(m)^2) - cosh(std_va(m)^2)) + ...
                         cos(va(m))^2*(sinh(2*std_va(m)^2) - sinh(std_va(m)^2))) + ...
                  std_vm(m)^2*exp(-2*std_va(m).^2)*...
                        (sin(va(m))^2*(2*cosh(2*std_va(m)^2) - cosh(std_va(m)^2)) + ...
                         cos(va(m))^2*(2*sinh(2*std_va(m)^2) - sinh(std_va(m)^2))));
                     
end

% standard deviation for current
NumMeasCurrent = length(i);

im = abs(i);
ia = unwrap(angle(i));

std_im = zeros(NumMeasCurrent,1);
std_ia = zeros(NumMeasCurrent,1);

std_i_re = zeros(NumMeasCurrent,1);
std_i_im = zeros(NumMeasCurrent,1);

for m = 1:NumMeasCurrent    
    
    if bus_label(m) == PCC
        rateI = 200/(baseA);
    else
        rateI = 40/(baseA);
    end
    
    % Linear Piece wise definition of current magnitude and phase meas std
    x_i = im(m)/rateI*100;
    
    if x_i < 0.2
        f_im = 15;
        f_ia = 300;
        
    else if x_i >= 0.2 && x_i < 1
            f_im = -10/0.8*x_i + 17.5;
            f_ia = -200/0.8*x_i + 350;
            
        else if  x_i >= 1 && x_i < 5
                f_im = -3.4/4*x_i + 5.85;
                f_ia = -72/4*x_i + 118;
                
            else if x_i >=5 && x_i < 20
                    f_im = -0.8/15*x_i + 1.8667;
                    f_ia = -13/15*x_i + 32.3333;
                    
                else if x_i >= 20 && x_i < 100
                        f_im = -0.3/80*x_i + 0.875;
                        f_ia = -6/80*x_i + 16.5;
                        
                    else if x_i >= 100
                            f_im = 0.5;
                            f_ia = 9;
                        end
                    end
                end
            end
        end
    end
    % std for abs and phase
    std_im(m) = (f_im/100)*im(m)/3;
    std_ia(m) = f_ia/3*1e-3;
    
    % std for real and imag
    std_i_re(m) = sqrt(...
                    im(m)^2*exp(-2*std_ia(m).^2)*...
                        (cos(ia(m))^2*(cosh(2*std_ia(m)^2) - cosh(std_ia(m)^2)) + ...
                         sin(ia(m))^2*(sinh(2*std_ia(m)^2) - sinh(std_ia(m)^2))) + ...
                    std_im(m)^2*exp(-2*std_ia(m).^2)*...
                        (cos(ia(m))^2*(2*cosh(2*std_ia(m)^2) - cosh(std_ia(m)^2)) + ...
                         sin(ia(m))^2*(2*sinh(2*std_ia(m)^2) - sinh(std_ia(m)^2))));
    
   std_i_im(m) = sqrt(...
                    im(m)^2*exp(-2*std_ia(m).^2)*...
                        (sin(ia(m))^2*(cosh(2*std_ia(m)^2) - cosh(std_ia(m)^2)) + ...
                         cos(ia(m))^2*(sinh(2*std_ia(m)^2) - sinh(std_ia(m)^2))) + ...
                    std_im(m)^2*exp(-2*std_ia(m).^2)*...
                        (sin(ia(m))^2*(2*cosh(2*std_ia(m)^2) - cosh(std_ia(m)^2)) + ...
                         cos(ia(m))^2*(2*sinh(2*std_ia(m)^2) - sinh(std_ia(m)^2))));
end

end

