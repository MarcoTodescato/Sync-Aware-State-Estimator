function [V, I, S, DeltaPhase, SE] = subsetEPFLdata(NumPMU, MeasIdx, SelectedPhase, PCC, DATA)
% Extract a desired subset of data from the preloaded set of EPFL data
%
% Args:
%   NumPMU: number of PMU
%   MeasIdx: indices of measurement to extract
%   SelectedPhase: phase to extract from the 3-phase data
%   PCC: index of the PCC node
%   DATA: preloaded data struct
%
% Returns:
%   V: complex voltage measurements
%   I: complex current measurements
%   S: complex current measurements
%   DeltaPhase: cumulative phase measurements
%   SE: state estimate data

NumIter = length(MeasIdx);

V = zeros(NumIter, NumPMU);
I = zeros(NumIter, NumPMU);
Frequency = zeros(NumIter, NumPMU);

% loop over PMU data
for id = 1:NumPMU
    
    % reference to data of 2014 Nov 17th 10-11h
    CurrentDataPMU = DATA.PMUs(id).Years.Months(11).Days(17).Hours(12);
    
    % convert in complex voltages
    V(:,id) = CurrentDataPMU.VM(MeasIdx, 1, SelectedPhase).*...
                   exp(1i.*(CurrentDataPMU.VA(MeasIdx, 1, SelectedPhase)));
    
    % convert in complex currents
    if id ~= PCC
        I(:,id) = -CurrentDataPMU.IM(MeasIdx,1,SelectedPhase).*...
                        exp(1i*(CurrentDataPMU.IA(MeasIdx, 1, SelectedPhase)));
    else
        I(:,id) = CurrentDataPMU.IM(MeasIdx,1,SelectedPhase).*...
                       exp(1i*(CurrentDataPMU.IA(MeasIdx, 1, SelectedPhase)));        
    end
    % frequency measurements
    Frequency(:,id) = CurrentDataPMU.Frequency(MeasIdx); % Hz
end

% cumulative phase
DeltaPhase = [0 ; cumsum(2*pi*(Frequency(1:end-1,1)./50 - 1), 1)];

% complex powers
S = V.*conj(I);

% 
CurrentDataSE = DATA.SE.Years.Months(11).Days(17).Hours(12);
SE.VM = CurrentDataSE.VM(MeasIdx, :, SelectedPhase);
SE.VA = CurrentDataSE.VA(MeasIdx, :, SelectedPhase);

end

