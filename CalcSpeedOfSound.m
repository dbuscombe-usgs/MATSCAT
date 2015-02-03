function [c] = CalcSpeedOfSound(T,S,D)

    if nargin ==0 
        fprintf('Temperature is required for Speed of Sound');
    elseif nargin==1
        S = 0;
        D = 0;
    elseif nargin==2
        D = 0;
    end
            
    c = 1448.96 + 4.591*T - (5.30e-2)*T^2 + (2.374e-4)*T^3 + 1.340*(S-35) + (1.630e-2)*D + (1.675e-7)*D^2 - (1.025e-2)*T*(S-35) - (7.139e-13)*T*D^3;
end 