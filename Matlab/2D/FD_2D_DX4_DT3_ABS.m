%% FD_2D_DX4_DT3_ABS.m 2-D acoustic Finite-Difference modelling
% GNU General Public License v3.0
%
% Author: Florian Wittkamp 2016
%
% Finite-Difference acoustic seismic wave simulation 
% Discretization of the first-order acoustic wave equation
%
% Temporal second-order accuracy O(DT^3)
% Spatial fourth-order accuracy  O(DX^4)
%
% Temporal discretization is based on the Adams-Basforth method
% Theory is available in:
% Bohlen, T., & Wittkamp, F. (2016). 
% Three-dimensional viscoelastic time-domain finite-difference 
% seismic modelling using the staggered Adams?Bashforth time integrator. 
% Geophysical Journal International, 204(3), 1781-1788.

%% Initialisation
disp(' ');
disp(['Starting ', mfilename ]);
close all, clear all;
seiscolor=load('colorbar.txt');

%% Input Parameter

% Discretization
c1=20;  % Number of grid points per dominant wavelength
c2=0.5; % CFL-Number
nx=200; ny=200; % Number of grid points
T=10;    % Total propagation time 

% Source Signal
f0= 5;      % Center frequency Ricker-wavelet
q0= 1;       % Maximum amplitude Ricker-Wavelet
xscr = 100;  yscr = 100; % Source position (in grid points)

% Receiver
xrec1=100; yrec1=120;  % Position Reciever 1 (in grid points)
xrec2=100; yrec2=130;  % Position Reciever 2 (in grid points)
xrec3=100; yrec3=140;  % Position Reciever 3 (in grid points)

% Velocity and density
modell_v=[3000*ones(ny,nx)]; % P-wave velocity in m/s
rho=[2200*ones(ny,nx)];        % Density in g/cm^3

%% Preparation 

% Init wavefields
vx=zeros(ny,nx);
vy=zeros(ny,nx);
p=zeros(ny,nx);

vx_x=zeros(ny,nx);
vy_y=zeros(ny,nx);
p_x=zeros(ny,nx);
p_y=zeros(ny,nx);


vx_x2=zeros(ny,nx);
vy_y2=zeros(ny,nx);
p_x2=zeros(ny,nx);
p_y2=zeros(ny,nx);

vx_x3=zeros(ny,nx);
vy_y3=zeros(ny,nx);
p_x3=zeros(ny,nx);
p_y3=zeros(ny,nx);

% Calculate first Lame-Paramter
l=rho.*modell_v.*modell_v;

cmin=min(modell_v(:)); % Lowest P-wave velocity
cmax=max(modell_v(:)); % Highest P-wave velocity
fmax=2*f0;             % Maximum frequency
dx=cmin/(fmax*c1);     % Spatial discretization (in m)
dy=dx;
dt=dx/(cmax)*c2;       % Temporal discretization (in s)
lampda_min=cmin/fmax;  % Smallest wavelength

% Create space and time vector
disp(['Model size: x:',num2str(dx*nx),' in m']);
disp(['Temporal discretization: ',num2str(dt),' s']);
disp(['Spatial discretization: ',num2str(dx),' m']);
disp(['Number of gridpoints per minimum wavelength: ',num2str(lampda_min/dx)]);

% Create space and time vector
x=0:dx:(dx*nx-dx); % Space vector
y=0:dx:(dx*ny-dx); % Space vector
t=0:dt:(T-dt);     % Time vector
nt=numel(t);       % Number of time steps

% Plotting model
figure
subplot(2,1,1)
imagesc(y,x,modell_v)
xlabel('Depth in m')
ylabel('VP in m/s')
set(gca,'FontSize',16)
title('Model')
subplot(2,1,2)
imagesc(y,x,rho)
xlabel('Depth in m')
ylabel('Density in g/cm^3')
set(gca,'FontSize',16)

% Source signal - Ricker-wavelet
tau=pi*f0*(t-1.5/f0);
q=q0*(1-2*tau.^2).*exp(-tau.^2);

% Plotting source signal
figure(2)
plot(t,q)
title('Source signal Ricker-Wavelet')
xlabel('Time in s')
ylabel('Amplitude')

% Init Seismograms
Seismogramm=zeros(3,nt); % 3 Seismogramme eingerichtet

% Calculation of some coefficients
c1=9/(8*dx);
c2=1/(24*dx);
c3=9/(8*dy);
c4=1/(24*dy);
c5=1/dx^3;
c6=1/dy^3;
c7=1/dx^2;
c8=1/dy^2;
c9=dt^3/24;

%% Time stepping
disp('Starting time stepping...');
tic;
for n=2:nt-1;
    
     % Update velocity
    for kx=5:nx-5;
        for ky=5:ny-5;
            
            % Calculating spatial derivative
                        p_x(ky,kx)=c1*(p(ky,kx+1)-p(ky,kx))...
                -c2*(p(ky,kx+2)-p(ky,kx-1));
            p_y(ky,kx)=c3*(p(ky+1,kx)-p(ky,kx))...
                -c4*(p(ky+2,kx)-p(ky-1,kx));
            
            % Update velocity
            vx(ky,kx)=vx(ky,kx)...
                -dt/rho(ky,kx)*1/24*(25*p_x(ky,kx)-2*p_x2(ky,kx)+p_x3(ky,kx));
            vy(ky,kx)=vy(ky,kx)...
                -dt/rho(ky,kx)*1/24*(25*p_y(ky,kx)-2*p_y2(ky,kx)+p_y3(ky,kx));
        end
    end
    
    vx_x3=vx_x2;
    vx_x2=vx_x;
    vy_y3=vy_y2;
    vy_y2=vy_y;
    
    % Inject source wavelet
    p(yscr,xscr)=p(yscr,xscr)+(q(n)+q(n-1))/2;
    
    % Update pressure
    for kx=5:nx-5;
        for ky=5:ny-5;
            
            % Calculating spatial derivative
            vx_x(ky,kx)=c1*(vx(ky,kx)-vx(ky,kx-1))...
                -c2*(vx(ky,kx+1)-vx(ky,kx-2));
            vy_y(ky,kx)=c3*(vy(ky,kx)-vy(ky-1,kx))...
                -c4*(vy(ky+1,kx)-vy(ky-2,kx));
            
             % Update velocity
            p(ky,kx)=p(ky,kx)...
                -l(ky,kx)*dt*1/24*(25*(vx_x(ky,kx)+vy_y(ky,kx))...
                -2*(vx_x2(ky,kx)+vy_y2(ky,kx))+(vx_x3(ky,kx)+vy_y3(ky,kx)));
        end
    end
    
    p_x3=p_x2;
    p_x2=p_x;
    p_y3=p_y2;
    p_y2=p_y;
    
    % Save seismograms 
    Seismogramm(1,n)=p(yrec1,xrec1);
    Seismogramm(2,n)=vy(yrec2,xrec2);
    Seismogramm(3,n)=p(yrec3,xrec3);
    
    % Plot Snapshot
    if(~mod(n,20))
        figure(3)
        imagesc(y,x,p)
        hold on;
        colormap([seiscolor])
        colorbar;
        set(gca,'CLim',[-0.05 0.05])
        set(3,'Position',[570 678 840 630])
        ylabel('Depth in m')
        xlabel('Horizontal in m')
        title(['Wavefield bei ', num2str(n*dt), 's'])
        rec1=plot(xrec1*dx-dx,yrec1*dy-dy);
        set(rec1,'MarkerSize',20,'Marker','.');
        rec2=plot(xrec2*dx-dx,yrec2*dy-dy);
        set(rec2,'MarkerSize',20,'Marker','.');
        rec3=plot(xrec3*dx-dx,yrec3*dy-dy);
        set(rec3,'MarkerSize',20,'Marker','.');
        scr=plot(xscr*dx-dx,yscr*dy-dy);
        set(scr,'MarkerSize',20,'Marker','.','Color',[1 0 0]);
        drawnow;
       
    end
end
time=toc;
N_t=numel(2:nt-1);
disp(['Finished with ', num2str(N_t),' time steps, in ',...
    num2str(time),' s']);
disp(['Time pro time step: ', num2str(time/N_t),' in s (Average)']);

%% Plot seismograms
figure(4)
set(4,'Position',[570 110 560 420])
subplot(3,1,1)
plot(t,Seismogramm(1,:));
xlabel('Time in s');
ylabel('Amplitude');
title(['Receiver at ', num2str(dx*yrec1), ' m']);
subplot(3,1,2)
plot(t,Seismogramm(2,:));
xlabel('Time in s');
ylabel('Amplitude');
title(['Receiver at ', num2str(dx*yrec2), ' m']);
subplot(3,1,3)
plot(t,Seismogramm(3,:));
xlabel('Time in s');
ylabel('Amplitude');
title(['Receiver at ', num2str(dx*yrec3), ' m']);

disp(' ');

% Save seismogram in matlab format
output=['Seismograms/',mfilename,'.mat'];
save(output,'Seismogramm','dt','T','time','dx');
