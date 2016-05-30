%% FD_1D_check_stability_func.m 1-D FD stability criterium check
% GNU General Public License v3.0
%
% Author: Florian Wittkamp 2016
%
% Check stability limit of FD-simulations
%
% Stability limit is calculated in terms of the CFL-number, which is
% defined as: CFL=v_(max)*DT/DX
% You get the maximum DT by DT=CFL*DX/v_(max)
%
% Theory:
% Fei, X., & Xiaohong, T. (2006).
% Stability and numerical dispersion analysis of a fourth-order
% accurate FDTD method. Antennas and Propagation,
% IEEE Transactions on, 54(9), 2525-2530.

%% Initialisation
function [bool]=FD_1D_check_stability_func(temporal_order,spatial_order,CFL)

%% Calculate spatial derivation impact
% Assumption: Spatial sampling is at the Nyquist condition
sum_fd_stencil=sum(abs(FD_taylor_coeff_func(spatial_order)));

switch temporal_order
    case 2
        % Temporal second-order (Leapfrog)
        theta=0:0.01:2*pi;
        f=-1i*(((exp(-1i*theta)).^2.5-(exp(-1i*theta)).^1.5));
        f2=(2*sum_fd_stencil*(1));
        f_ges=abs(f./f2);
        if(f_ges<CFL)
            bool=0;
        else
            bool=1;
        end
        
    case 3
        % Temporal third-order (Adams-Bashforth method)
        theta=0:0.01:2*pi;
        f=-1i*(((exp(-1i*theta)).^2.5-(exp(-1i*theta)).^1.5)*1);
        f2=(2*sum_fd_stencil.*(25/24.*(exp(-1i*theta)).^2-1/12.*(exp(-1i*theta)).^1+1/24));
        f_ges=abs(f./f2);
        if(f_ges<CFL)
            bool=0;
        else
            bool=1;
        end
        
    case 4
        % Temporal fourth-order (Adams-Bashforth method)
        theta=0:0.01:2*pi;
        f=-1i*(((exp(-1i*theta)).^3.5-(exp(-1i*theta)).^2.5)*1);
        f2=(2*sum_fd_stencil.*(13/12.*(exp(-1i*theta)).^3-5/24.*(exp(-1i*theta)).^2+1/6.*(exp(-1i*theta)).^1-1/24));
        f_ges=abs(f./f2);
        if(f_ges<CFL)
            bool=0;
        else
            bool=1;
        end
        
    otherwise
        error('Temporal-order not supported');
end
end