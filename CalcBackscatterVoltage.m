function [V,r] =  CalcBackscatterVoltage(Model,Freq,At,Kt,r,M,As,deltaAs,T,S,Depth)
%--------------------------------------------------------------------------
% 
% The purpose of this function is to generate a theoretical backscatter
% voltage given a set of environmental parameters. This function assumes
% that the bins are a consistent size which is typical of such a system.
%
%   Inputs
%   Model = 'GLASS', 'SAND'
%   Freq in Hz
%   At is equivalent Transducer Radius
%   Kt is the system constant
%   r[N] is the range vector (m) where N is the number of bins 
%   M[N] is the mass concentration with range (kgm-1)
%   As{N] is the sediment radius in m with range  
%   Delta As in m
%   Density is the Sediment density
%   T  is temperature in degC
%   S  is Salinity in ppm
%   Depth is the instrument depth in m
%--------------------------------------------------------------------------
Depth = 0;
BinSize = r(2) - r(1);
NumBins = size(r,2);
%Speed of Sound
c = CalcSpeedOfSound(T,S,Depth);
r = r * c/1500;     % Correct for speed of sound
%Water Attenuation in nepers/m assume constant temperature profile
AttenWater = CalcWaterAttenuation(Freq,T,S,Depth);

%If Mass or Sediment Size is not an array then assume same for all Bins
if size(M,1) == 1 
    M(1:NumBins,1) = M;
end
if size(As,1) == 1 
    As(1:NumBins,1) = As;
end

% Deal with the sediment form function and scattering cross section
% assumes fixed size sediment distribution
[ff XX Density] = CalcFormFunction(Model,As,deltaAs,Freq, c);

%Deal with sediment attenuation in each bin
AttenSed = (M .* ((4*XX)./(3*Density*As)) * BinSize); 
AttenSed(1) = AttenSed(1)*r(1)/BinSize;        %Deal with offset first bin
AttenSed = cumsum(AttenSed);                   %Now Accumulate so range dependent

%Range Correction
Rn = (pi * (At.^2) .* Freq)/ c;
z = r./Rn;
Yeta = (1+1.35*z+((2.5*z).^3.2))./(1.35*z+((2.5*z).^3.2));

%Calculate Sediment Constant Ks
Ks = ff./sqrt(As*Density);


V = ((Kt * Ks) .* M.^0.5 .* exp(-2*r*AttenWater) .* exp(-2*AttenSed)) ./ (r.*Yeta);

end



