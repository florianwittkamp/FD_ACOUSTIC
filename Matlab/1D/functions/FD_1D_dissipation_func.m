%% FD_1D_dissipation_func.m Calculate dissipation for 1D acoustic FD-schemes
% GNU General Public License v3.0
%
% Author: Florian Wittkamp 2016
%
% Calculating the numerical dissipation as remaning amplitude after one
% time step.
%
% Supports only the second-order leapfrog scheme and the ABS third and
% fourth-order schemes.

function [return_value]=FD_1D_dissipation_func(temporal_order,spatial_order,CFL,KH)

%% Calculate spatial derivation impact
coeff_spatial=FD_taylor_coeff_func(spatial_order);
N=numel(coeff_spatial);
spatial_influence=0;
for n=1:N;
    spatial_influence=spatial_influence+str2num(rats(coeff_spatial(n)))*sin(KH/2*(2*n-1));
end

%% Calculate dissipation for a single time step
switch temporal_order;
    case 2;
        % The second-order leapfrog scheme does not suffer from dissipation
        dissipation=1;
    case 3;
        if(FD_1D_check_stability_func(temporal_order,spatial_order,CFL)==0)
            dissipation=NaN;
        else
            syms t;
            dissipation=vpasolve(t^(5/2)-t^(3/2)-2*1i*CFL*spatial_influence*(25/24*t^2-1/12*t+1/24)==0,t,1);
            dissipation=abs(dissipation);
        end
    case 4;
        if(FD_1D_check_stability_func(temporal_order,spatial_order,CFL)==0)
            dissipation=NaN;
        else
            syms t;
            dissipation=vpasolve(t^(7/2)-t^(5/2)-2*1i*CFL*spatial_influence*(13/12*t^3-5/24*t^2+1/6*t-1/24)==0,t,1);
            dissipation=abs(dissipation);
        end
    otherwise
        error('Temporal-order not implemented')
end
return_value=dissipation;
end
