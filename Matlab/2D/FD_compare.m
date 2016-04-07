%% FD_compare.m 2-D acoustic Finite-Difference modelling
% GNU General Public License v3.0
%
% Author: Florian Wittkamp 2016
%
% Finite-Difference acoustic seismic wave simulation 
% Discretization of the first-order acoustic wave equation
%
% Compare seismograms, which are calculated with different 
% spatial and temporal accuracy.
%
% This script will run each FD-code and read in the seismograms to
% compare them. 

%% Initialisation
disp(' ');
disp(['Starting ', mfilename ]);
close all, clear all;

%% Calculate Seismograms:

% Call each FD script, each will save the seismogram to Seismograms/
% This will need some time...
FD_2D_DX4_DT2
FD_2D_DX4_DT3_ABS

% Clear memory and close all figures
close all; clear all;

%% Read Seismograms
FD_4_2=load('Seismograms/FD_2D_DX4_DT2.mat');
FD_4_3_ABS=load('Seismograms/FD_2D_DX4_DT3_ABS.mat');

% Calculate time vector
t=0:FD_4_2.dt:(FD_4_2.T-FD_4_2.dt);

%% Plotting
figure
plot(t,FD_4_2.Seismogramm(1,:))
hold on
plot(t,FD_4_3_ABS.Seismogramm(1,:))
set(gca,'FontSize',16)
xlabel('Time in s')
ylabel('Amplitude')
title('Seismogram 1')
legend('FD 4 2','FD 4 3 ABS')

figure
plot(t,FD_4_2.Seismogramm(2,:))
hold on
plot(t,FD_4_3_ABS.Seismogramm(2,:))
set(gca,'FontSize',16)
xlabel('Time in s')
ylabel('Amplitude')
title('Seismogram 2')
legend('FD 4 2','FD 4 3 ABS')

figure
plot(t,FD_4_2.Seismogramm(3,:))
hold on
plot(t,FD_4_3_ABS.Seismogramm(3,:))
set(gca,'FontSize',16)
xlabel('Time in s')
ylabel('Amplitude')
title('Seismogram 3')
legend('FD 4 2','FD 4 3 ABS')


disp(' ');