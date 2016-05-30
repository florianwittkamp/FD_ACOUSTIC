%% FD_taylor_coeff_func.m calculate Taylor coefficient
% GNU General Public License v3.0
%
% Author: Florian Wittkamp 2016
%
% Calculate the Taylor coefficient in an arbitrary order for the
% second-order derivative.
%
% Usage:
% Lets say you want to calculate the 4th order accurate second-order
% FD-stencil. Then you have to set order=4 and the result would be used
% as follow:
% p_x = 1/DH * ( coeff(1) * (p(x+1)-p(x)) + coeff(2) * (p(x+2)-p(x-1)) )
% where p_x is the derivative.

function coeff=FD_taylor_coeff_func(spatial_order)

%% Check some things:
if(mod(spatial_order,2)~=0)
    disp(num2str(spatial_order))
    error('Spatial-order has to be an integer multiple of 2!')
    % There will be only derivation with an order that is an integer multiple of 2
end

if(spatial_order<2)
    error('Order has to be at least 2!')
end

%% Calculation

if(spatial_order==2)
    % Return result
    coeff=1;
end

if(spatial_order>2)
    c=[1 zeros(1,spatial_order/2-1)]';
    M=zeros(2,spatial_order/2);
    % Condition 1: \sum^{N/2}_{k=1} b_k(2k-1)=1
    for n=1:spatial_order/2;
        M(1,n)=(2*n-1);
    end
    % Condition 2:  \sum^{N/2}_{k=1} b_k(2k-1)^(2j-1)=0; j=2,3...N/2
    for j=2:spatial_order/2;
        for n=1:spatial_order/2;
            M(j,n)=(2*n-1)^(2*j-1);
        end
    end
    % Return result
    coeff=(M\c)';
end

end