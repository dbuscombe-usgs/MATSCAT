function [atten] = CalcWaterAttenuation(Freq, T, S, D)
%-------------------------------------------------------------------------
%   Function to calculate the water attenuation for a specific frequency
%   using the ???? Model. 
%
%   Parameters
%   Freq = Acoustic Frequency in Hz
%   T = Temperature in degC
%   S = Salinity in ppm (default = 0)
%   D = Depth in metres (default = 0)
%   returns attenuation in Nepers/m
%-------------------------------------------------------------------------
%   Version History
%   Version 1.0 Only uses temperature as main contributor above 1MHz
%-------------------------------------------------------------------------

   if nargin < 2 
        fprintf('Requires Frequency and Temperature');
   elseif nargin==2
        S = 0;
        D = 0;
   elseif nargin==3
        D = 0;
   end
           

   % atten = ((2.1e-10)*( T-38)^2+(1.3e-7))*(((Freq/1000).^2) / 8.686)	%In Nepers
  
  % Equation below is Attenuation for fresh water calculated using Fisher and Simmons
   F = Freq/1E6;
   atten = ((55.9 - 2.37*T + 4.77E-2*T^2 - 3.84E-4*T^3)*10^-3*F^2);
end
