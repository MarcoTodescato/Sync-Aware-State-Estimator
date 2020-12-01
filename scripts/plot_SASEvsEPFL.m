% Script to plot result of SASE vs Algorithm presented in
%   M. Pignati, et al. Real-Time State Estimation of the EPFL-Campus Medium-Voltage Grid by Using PMUs. 
%   The Sixth Conference on Innovative Smart Grid Technologies (ISGT2015), 2015.
% on the EPFL campus 5 node testbed

close all
clc

% To replicate paper results load corresponding data
%load('./data/publications/Energies_SASEvsEPFL.mat');

% reload temporal info for plot
T = 1/PMURate:1/PMURate:SimTime;
MeasIdxVisualization = (FromMin*60 + FromSec)*PMURate + (1:PMURate*SimTime)';
Idx = ismember(MeasIdx,MeasIdxVisualization);

LineWidth = 1.5;
FontSizeLegend = 12;
FontSizeLabels = 14;
FontSizeAxis = 11;


%--------%
% powers %
%--------%
% measured
switch measurement_type

    case 'voltage+power'
        Rpowers = R(2*length(LabelPMUMeas)+1:end,2*length(LabelPMUMeas)+1:end,:);

    otherwise
        Rpowers = zeros(2*length(LabelPMUMeas));
end

for j = 3 % 1:length(LabelPMUMeas)

    id = LabelPMUMeas(j);

    figure
    subplot(2,1,1)
    hold on
    switch measurement_type
        case 'voltage+power'
            fill([T,flip(T)],[pSASE(Idx,id)+3*sqrt(squeeze(PpSASE(id,id,:)))+3*sqrt(squeeze(Rpowers(2*j-1,2*j-1,:)));...
              flip(pSASE(Idx,id)-3*sqrt(squeeze(PpSASE(id,id,:)))-3*sqrt(squeeze(Rpowers(2*j-1,2*j-1,:))))],...
             '-','FaceColor',rgb('Cyan'),'FaceAlpha',0.1,'EdgeColor','none'); 
        otherwise
    end
    fill([T,flip(T)],[pSASE(Idx,id)+3*sqrt(squeeze(PpSASE(id,id,:)));...
         flip(pSASE(Idx,id)-3*sqrt(squeeze(PpSASE(id,id,:))))],...
         '-','FaceColor',rgb('Red'),'FaceAlpha',0.1,'EdgeColor','none');
    plot(T,real(S(Idx,id)),'k-','LineWidth',LineWidth)
    plot(T,pEPFL(Idx,id),'--','Color',rgb('DarkGreen'),'LineWidth',LineWidth)
    plot(T,pSASE(Idx,id),'-','Color',rgb('Red'),'LineWidth',LineWidth)
    switch measurement_type
        case 'voltage+power'
            hl = legend('Confidence + Meas. Cov','Confidence','Measurements','EPFL','SASE');
        otherwise
            hl = legend('Confidence','Measurements','Alg. [12]','SASE');
    end
    set(hl,'FontSize',FontSizeLegend,'Location','SouthWest')
    set(gca,'FontSize',FontSizeAxis);
    hy1 = ylabel('P');
    hold off
    grid on
    xlim([2 SimTime])
    subplot(2,1,2)
    hold on
    switch measurement_type
        case 'voltage+power'
            fill([T,flip(T)],[qSASE(Idx,id)+3*sqrt(squeeze(PqSASE(id,id,:)))+3*sqrt(squeeze(Rpowers(2*j,2*j,:)));...
                flip(qSASE(Idx,id)-3*sqrt(squeeze(PqSASE(id,id,:)))-3*sqrt(squeeze(Rpowers(2*j,2*j,:))))],...
                '-','FaceColor',rgb('Cyan'),'FaceAlpha',0.1,'EdgeColor','none');
        otherwise
    end
    fill([T,flip(T)],[qSASE(Idx,id)+3*sqrt(squeeze(PqSASE(id,id,:)));...
        flip(qSASE(Idx,id)-3*sqrt(squeeze(PqSASE(id,id,:))))],...
       '-','FaceColor',rgb('Red'),'FaceAlpha',0.1,'EdgeColor','none');
    plot(T,imag(S(Idx,id)),'k-','LineWidth',LineWidth)
    plot(T,qEPFL(Idx,id),'--','Color',rgb('DarkGreen'),'LineWidth',LineWidth)
    plot(T,qSASE(Idx,id),'-','Color',rgb('Red'),'LineWidth',LineWidth)
    hy2 = ylabel('Q');
    hx = xlabel('Time[s]');
    hold off
    grid on
    set(hy1,'FontSize',FontSizeLabels,'Interpreter','latex')
    set(hy2,'FontSize',FontSizeLabels,'Interpreter','latex')
    set(hx,'FontSize',FontSizeLabels,'Interpreter','latex')
    set(gca,'FontSize',FontSizeAxis);
    xlim([2 SimTime])

end

% predicted
LabelPMUPred = setdiff(LabelPMU,LabelPMUMeas);
for id = 5; % LabelPMUPred

    figure
    subplot(2,1,1)
    hold on
    fill([T,flip(T)],[pSASE(Idx,id)+3*sqrt(squeeze(PpSASE(id,id,:)));...
        flip(pSASE(Idx,id)-3*sqrt(squeeze(PpSASE(id,id,:))))],...
        '-','FaceColor',rgb('Red'),'FaceAlpha',0.1,'EdgeColor','none');
    plot(T,real(S(Idx,id)),'k-','LineWidth',LineWidth)
    %plot(T,pEPFL(Idx,id),'--','Color',rgb('DarkGreen'),'LineWidth',LineWidth)
    plot(T,pSASE(Idx,id),'-','Color',rgb('Red'),'LineWidth',LineWidth)
    %hl = legend('Confidence','Measurements','EPFL','SASE');
    hl = legend('Confidence','Measurements','SASE');
    set(hl,'FontSize',FontSizeLegend,'Location','NorthWest')
    set(gca,'FontSize',FontSizeAxis);
    hy1 = ylabel('P');
    hold off
    grid on
    xlim([2 SimTime])
    subplot(2,1,2)
    hold on
    fill([T,flip(T)],[qSASE(Idx,id)+3*sqrt(squeeze(PqSASE(id,id,:)));...
        flip(qSASE(Idx,id)-3*sqrt(squeeze(PqSASE(id,id,:))))],...
        '-','FaceColor',rgb('Red'),'FaceAlpha',0.1,'EdgeColor','none');
    plot(T,imag(S(Idx,id)),'k-','LineWidth',LineWidth)
    %plot(T,qEPFL(Idx,id),'--','Color',rgb('DarkGreen'),'LineWidth',LineWidth)
    plot(T,qSASE(Idx,id),'-','Color',rgb('Red'),'LineWidth',LineWidth)
    hy2 = ylabel('Q');
    hx = xlabel('Time[s]');
    hold off
    grid on
    set(hy1,'FontSize',FontSizeLabels,'Interpreter','latex')
    set(hy2,'FontSize',FontSizeLabels,'Interpreter','latex')
    set(hx,'FontSize',FontSizeLabels,'Interpreter','latex')
    set(gca,'FontSize',FontSizeAxis);
    xlim([2 SimTime])

end

