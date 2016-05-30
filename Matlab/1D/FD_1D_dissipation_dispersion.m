%% Plot dissipation and dispersion for a 1D acoustic FD-scheme
% GNU General Public License v3.0
%
% Author: Florian Wittkamp 2016
%
% Plotting the numerical dissipation and numerical dispersion for
% the 1D acoustic FD-schemes.
% This script supports the second-order leapfrog scheme as well as
% the third and fourth-order Adams-Basforth schemes.

%% Initialisation
disp(' ');
disp(['Starting ', mfilename ]);
clearvars; close all;
addpath functions

%% Input Parameter

Courant_Max=0.8; % Maximum Courant (CFL) number
Courant_Min=0.01; % Smallest Courant (CFL) number
Courant_Delta=0.05; % Sampling for the Courant (CFL) number
KH=2*pi/15; % Spatial sampling, 2*pi/(Grid points per wavelength)
Spatial_order=8; % Spatial-order of the FD-stencil
NT=2000; % Number of time steps for the dissipation calculation

%% Calculating the numerical dissipation and dispersion
CFL=Courant_Min:Courant_Delta:Courant_Max;
for order=2:1:4;
    n=1;
    disp(['Calculating dissipation/dispersion for temporal-order: ',num2str(order)]);
     for c=CFL;
        Dissipation(order-1,n)=FD_1D_dissipation_func(order,Spatial_order,c,KH).^NT;
        Dispersion(order-1,n)=FD_1D_dispersion_func(order,Spatial_order,c,KH);
        n=n+1;
    end
end
disp('Calculation finished');

%% Plotting of the dispersion and dissipation
subplot(1,2,1)
p=plot(CFL,(Dissipation),'LineWidth',2);
legend('M=2','M=3','M=4')
title('Dissipation')
xlabel('CFL-Number')
ylabel(sprintf(['Amplitude after \n',num2str(NT),'time steps']))
set(gca,'FontSize',16)
set(p(1),'color',[0.8500 0.3250 0.0980])
set(p(2),'color',[0.9290 0.6940 0.1250])
set(p(3),'color',[0.4940 0.1840 0.5560])

subplot(1,2,2)
p=plot(CFL,Dispersion-1,'LineWidth',2);
legend('M=2','M=3','M=4')
title('Dispersion')
xlabel('CFL-Number')
ylabel('c_{fd}/c-1')
set(gca,'FontSize',16)
hold on
one=zeros(1,numel(CFL));
h2=plot(CFL,one,'-.','color','black');
set(p(1),'color',[0.8500 0.3250 0.0980])
set(p(2),'color',[0.9290 0.6940 0.1250])
set(p(3),'color',[0.4940 0.1840 0.5560])