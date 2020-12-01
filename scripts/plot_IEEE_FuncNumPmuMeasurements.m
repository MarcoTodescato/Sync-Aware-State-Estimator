% Scripts to plot the performance of SASE tested on IEEE synthetic
% loadcases as function of the number measurements per PMU (and a fixed
% number of PMU deployed)
%
% NOTE: if plotting your own data, axis limits might need to be turned off
%
% NOTE: when running the script standalone, make sure to load data  saved 
% by running main_SASEonIEEE.m with the correctly set execution mode:
%   exe_mode = mode_1

close all
clc

addpath(genpath('./'));
 
% To replicate paper results load corresponding data
%load('./data/publications/Energies_IEEE123_FuncNumPmuMeasurements.mat');
%load('./data/publications/CDC_CASE15_FuncNumPmuMeasurements.mat');

LineWidth = 2;
FontSizeAxis = 14;
FontSizeLegend = 12;
FontSizeLabel = 15;

Colors = {'Green','DodgerBlue','Crimson','DarkOrange','DodgerBlue','DarkGreen'};

% ARMSE voltages
figure
hold on
plot(P.NUM_MEAS,log10(PFK.ERR.V),'Color',rgb(Colors{1}),'LineStyle','-','LineWidth',LineWidth);
plot(P.NUM_MEAS,log10(PFK.TR.V),'Color',rgb(Colors{1}),'LineStyle','--','LineWidth',LineWidth);
plot(P.NUM_MEAS,log10(SE.ERR.V),'Color',rgb(Colors{2}),'LineStyle','-','LineWidth',LineWidth)
plot(P.NUM_MEAS,log10(SE.TR.V),'Color',rgb(Colors{2}),'LineStyle','--','LineWidth',LineWidth)
plot(P.NUM_MEAS,log10(ABN.ERR.V),'Color',rgb(Colors{3}),'LineStyle','-','LineWidth',LineWidth)
plot(P.NUM_MEAS,log10(ABN.TR.V),'Color',rgb(Colors{3}),'LineStyle','--','LineWidth',LineWidth)
hold off
grid on
hl = legend('GT - Empirical','GT - Theoretical','SASE - Empirical','SASE - Theoretical',...
            'BLSE - Empirical','BLSE - Theoretical');
hy = ylabel('${\rm log}_{10}({\rm ARMSE}(\widehat{\mathbf{u}},\mathbf{u},\cdot)$)');
hx = xlabel('M');
xlim([1,P.NUM_MEAS(end)]);
pbaspect([1 .5 1]);
set(gca,'FontSize',FontSizeAxis);
set(hl,'FontSize',FontSizeLegend,'Location','East');
set(hy,'FontSize',FontSizeLabel,'Interpreter','Latex');
set(hx,'FontSize',FontSizeLabel,'Interpreter','Latex');


% ARMSE delay parameters
figure
hold on 
plot(P.NUM_MEAS,log10(SE.ERR.SKEW),'Color',rgb(Colors{4}),'LineStyle','-','LineWidth',LineWidth);
plot(P.NUM_MEAS,log10(SE.TR.SKEW),'Color',rgb(Colors{6}),'LineStyle','--','LineWidth',LineWidth);
plot(P.NUM_MEAS,log10(Std(3,:)),'Color',rgb(Colors{4}),'LineStyle',':','LineWidth',LineWidth);
plot(P.NUM_MEAS,log10(SE.ERR.OFF),'Color',rgb(Colors{5}),'LineStyle','-','LineWidth',LineWidth);
plot(P.NUM_MEAS,log10(SE.TR.OFF),'Color',rgb(Colors{3}),'LineStyle','--','LineWidth',LineWidth);
plot(P.NUM_MEAS,log10(Std(1,:)),'Color',rgb(Colors{5}),'LineStyle',':','LineWidth',LineWidth);
hold off
grid on
hl = legend('Skew - Empirical','Skew - Theoretical','Skew - Theoretical (2 nodes case)',...
            'Offset - Empirical','Offset - Theoretical','Offset - Theoretical (2 nodes case)');
hy = ylabel({'${\rm log}_{10}({\rm ARMSE}(\widehat{\alpha},\alpha,\cdot))$';...
    '${\rm log}_{10}({\rm ARMSE}(\widehat{\beta},\beta,\cdot))$'});
hx = xlabel('M');
pbaspect([1 .5 1]);
xlim([1,P.NUM_MEAS(end)]);
set(gca,'FontSize',FontSizeAxis);
set(hl,'FontSize',FontSizeLegend);
set(hy,'FontSize',FontSizeLabel,'Interpreter','Latex');
set(hx,'FontSize',FontSizeLabel,'Interpreter','Latex');

