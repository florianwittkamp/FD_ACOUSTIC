%% FD_1D_dispersion_func.m Calculate dispersion for 1D acoustic FD-schemes
% GNU General Public License v3.0
%
% Author: Florian Wittkamp 2016
%
% Calculating the numerical dispersion as ratio between numerical
% wave velocity and model velocity.
%
% Supports only the second-order leapfrog scheme and the ABS third and
% fourth-order schemes.

function [return_value]=FD_1D_dispersion_func(temporal_order,spatial_order,CFL,KH)

%% Calculate spatial derivation impact
coeff_spatial=FD_taylor_coeff_func(spatial_order);
N=numel(coeff_spatial);
spatial_influence=0;
for n=1:N;
    spatial_influence=spatial_influence+str2num(rats(coeff_spatial(n)))*sin(KH/2*(2*n-1));
end

%% Calculate dispersion
switch temporal_order;
    case 2;
        if(FD_1D_check_stability_func(temporal_order,spatial_order,CFL)==0)
            dispersion=NaN;
        else
            dispersion=asin(CFL*spatial_influence).*2./(KH*CFL);
        end
    case 3;
        if(FD_1D_check_stability_func(temporal_order,spatial_order,CFL)==0)
            dispersion=NaN;
        else
            syms c;
            dispersion=vpasolve(sin(0.5*c*CFL*KH)/(25/24-1/12*cos(1*c*CFL*KH)+1/24*cos(2*c*CFL*KH))-CFL*spatial_influence==0,c,1);
        end
    case 4;
        if(FD_1D_check_stability_func(temporal_order,spatial_order,CFL)==0)
            dispersion=NaN;
        else
            syms c;
            dispersion=vpasolve(sin(0.5*c*CFL*KH)/(13/12-5/24*cos(1*c*CFL*KH)+1/6*cos(2*c*CFL*KH)-1/24*cos(3*c*CFL*KH))-CFL*spatial_influence==0,c,1);
        end
    otherwise
        error('Temporal-order not implemented')
end
return_value=dispersion;
end