% Scripts to plot the performance of SASE tested on IEEE synthetic
% loadcases as function of the number of PMU deployed (and a fixed number
% of measurements per PMU)
%
% NOTE: plotting your own data, axis limits might need to be turned off
%
% NOTE: when running the script standalone, make sure to load data  saved 
% by running main_SASEonIEEE.m with the correctly set execution mode:
%   exe_mode = mode_2

close all
clc

addpath(genpath('./'));
 
% To replicate paper results load corresponding data
%load('./data/publications/Energies_IEEE123_FuncNumPmuDeployed.mat');
%load('./data/publications/CDC_CASE15_FuncNumPmuDeployed.mat');

LineWidth = 2;
FontSizeAxis = 14;
FontSizeLegend = 12;
FontSizeLabel = 15;

Colors = {'Green','DodgerBlue','Crimson','DarkOrange','DodgerBlue','DarkGreen'};

% ARMSE voltages
figure
hold on
plot(P.NUM_PMU,log10(PFK.ERR.V),'Color',rgb(Colors{1}),'LineStyle','-','LineWidth',LineWidth);
plot(P.NUM_PMU,log10(PFK.TR.V),'Color',rgb(Colors{1}),'LineStyle','--','LineWidth',LineWidth);
plot(P.NUM_PMU,log10(SE.ERR.V),'Color',rgb(Colors{2}),'LineStyle','-','LineWidth',LineWidth);
plot(P.NUM_PMU,log10(SE.TR.V),'Color',rgb(Colors{2}),'LineStyle','--','LineWidth',LineWidth);
plot(P.NUM_PMU,log10(ABN.ERR.V),'Color',rgb(Colors{3}),'LineStyle','-','LineWidth',LineWidth);
plot(P.NUM_PMU,log10(ABN.TR.V),'Color',rgb(Colors{3}),'LineStyle','--','LineWidth',LineWidth);
grid on
hold off
ylim([-4.3 -2])
hl = legend('GT - Empirical','GT - Theoretical','SASE - Empirical','SASE Theoretical',...
                'BLSE - Empirical','BLSE Theoretical');
hy = ylabel('${\rm log}_{10}({\rm ARMSE}(\widehat{\mathbf{u}},\mathbf{u},M))$');
hx = xlabel('Number of PMUs');
pbaspect([1 .5 1]);
set(gca,'FontSize',FontSizeAxis);
set(hl,'FontSize',FontSizeLegend,'Location','East');
set(hy,'FontSize',FontSizeLabel,'Interpreter','Latex');
set(hx,'FontSize',FontSizeLabel,'Interpreter','Latex');

% ARMSE delay parameters
figure
hold on
plot(P.NUM_PMU(2:end),log10(SE.ERR.SKEW(2:end,:)),'Color',rgb(Colors{4}),'LineStyle','-','LineWidth',LineWidth);
plot(P.NUM_PMU(2:end),log10(SE.TR.SKEW(2:end,:)),'Color',rgb(Colors{4}),'LineStyle','--','LineWidth',LineWidth);
plot(P.NUM_PMU(2:end),log10(SE.ERR.OFF(2:end,:)),'Color',rgb(Colors{5}),'LineStyle','-','LineWidth',LineWidth);
plot(P.NUM_PMU(2:end),log10(SE.TR.OFF(2:end,:)),'Color',rgb(Colors{5}),'LineStyle','--','LineWidth',LineWidth);
hold off
grid on
xlim([1,P.MaxNumPmu])
hl = legend('Skew - Empirical','Skew - Theoretical','Offset - Empirical','Offset - Theoretical');
hy = ylabel({'${\rm log}_{10}({\rm ARMSE}(\widehat{\alpha},\alpha,M))$';'${\rm log}_{10}({\rm ARMSE}(\widehat{\beta},\beta,M))$'});
hx = xlabel('Number of PMUs');
pbaspect([1 .5 1]);
set(gca,'FontSize',FontSizeAxis);
set(hl,'FontSize',FontSizeLegend,'Location','East');
set(hy,'FontSize',FontSizeLabel,'Interpreter','Latex');
set(hx,'FontSize',FontSizeLabel,'Interpreter','Latex');

