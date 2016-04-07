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

function coeff=FD_coeff_func(order);
%% Check some things:
if(mod(order,2)~=0)
    error('Order has to be a integer multiple of 2!')
    % There will by only derivative with a order which is 
    % a integer multiple of 2
end
if(order<=2)
    error('Order has to be at least 4!')
end
%% Calculation
c=[1 zeros(1,order/2-1)]';
M=zeros(2,order/2);
% Condition 1: \sum^{N/2}_{k=1} b_k(2k-1)=1
for n=1:order/2;
    M(1,n)=(2*n-1);
end
% Condition 2:  \sum^{N/2}_{k=1} b_k(2k-1)^(2j-1)=0; j=2,3...N/2
for j=2:order/2;
    for n=1:order/2;
        M(j,n)=(2*n-1)^(2*j-1);
    end
end
%% Return result
coeff=(inv(M)*c)';
end