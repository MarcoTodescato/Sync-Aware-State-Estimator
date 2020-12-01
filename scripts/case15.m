function mpc = case15
%--------------------------------------------------------------------------
% MATPOWER load case of the 15 nodes distribution feeder from 
%
% Simple and efficient method for load flow solution of radial distribution
% networks.
% D. Das, D.P. Kothari, A. Kalam.
% Electrical Power and Energy Systems 1995
% URL: http://www.sciencedirect.com/science/article/pii/0142061595000500
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% loadcase version
%--------------------------------------------------------------------------
mpc.version = '2';

%--------------------------------------------------------------------------
% base MVA
%--------------------------------------------------------------------------
mpc.baseMVA = 100;

%--------------------------------------------------------------------------
% bus
%--------------------------------------------------------------------------
mpc.bus = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; ...
           3 1 1 1 1 1 1 1 1 1 1 1 1 1 1; ...
           [0 48.1851 76.4842 152.9684 48.1851 152.9684 152.9684 76.4842 ...
           76.4842 48.1851 152.9684 76.4842 48.1851 76.4842 152.9684]./1000; ...
           [0 40.5857 64.4218 128.8435 40.5857 128.8435 128.8435 64.4218 64.4218 ...
           40.5857 128.8435 64.4218  40.5857 64.4218 128.8435]./1000; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
           1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; ...
           1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; ...
           0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
           1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; ...
           1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06
           0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94 0.94]';

%--------------------------------------------------------------------------
% generator
%--------------------------------------------------------------------------
mpc.gen = [1; ...
           1340/1000; ...
           1128.7/1000; ...
           2*1128.7/1000; ...
           -2*1128.7/1000; ...
           1; ...
           100; ...
           1; ...
           2*1340/1000; ...
           0; ...
           0; ...
           0; ...
           0; ...
           0; ...
           0; ...
           0; ...
           0; ...
           0; ...
           0; ...
           0; ...
           0]';

%--------------------------------------------------------------------------
% branch
%--------------------------------------------------------------------------
mpc.branch = [1 2 3 4 2 9 2 6 6 3 11 12 4 4; ...
              2 3 4 5 9 10 6 7 8 11 12 13 14 15; ...
              [1.35309 1.17024 0.84111 1.52348 2.01317 1.68671 2.55727 ...
              1.08820 1.25143 1.79553 2.44845 2.01317 2.23081 1.19702]./1.2; ...
              [1.32349 1.14464 0.82271 1.02760 1.35790 1.13770 1.72490 ...
              0.73400 0.84410 1.21110 1.65150 1.35790 1.50470 0.80740]./1.2; ...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
              0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
              1 1 1 1 1 1 1 1 1 1 1 1 1 1; ...
              -360 -360 -360 -360 -360 -360 -360 -360 -360 -360 -360 -360 -360 -360; ...
              360 360 360 360 360 360 360 360 360 360 360 360 360 360]';

