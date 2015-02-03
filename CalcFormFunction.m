function [ffa xxa Density] = CalcFormFunction(Model, As,DevAs, f, c)
%--------------------------------------------------------------------------
%   [ffa xxa Density] = CalcFormFunction(Model, As,DevAs, f, c)
% 
%   Calculates the required Form Function and scattering cross section for
%   a given frequency. The module takes into account the size distribution. 
%   The module will cope with an array of sediment size and frequency
%   It uses the CalcSievedFormFunction routine to determine the
%   form function and sets the min and max limits to 3x the standard
%   deviation. 
%--------------------------------------------------------------------------
%   Passed
%   Model = GLASS or SAND
%   As = Mean Sediment Radius (m) (Can be array)
%   DevAs = Deviation of Sediment Size (sigma = DevAs*As)
%   f = acoustic frequency in Hz [F1 F2 F3 ..]
%   c = speed of sound in water in m/s
%   Returns
%   ffa[NumSizes,NumFreq] = form function organised as ff[numSizes,numFreq]
%   xxa[NumSizes,NumFreq] = scattering cross-section as xx[numSizes,numFreq]
%   Density = Denisty of typical material
%-------------------------------------------------------------------------
%   Version
%   1.0     Initial Version
%   2.0     Rewritten to use the CalcSievedFormFunction 
%   2.1     Corrected syntax errors
%------------------------------------------------------------------------- 

%Calculate the MinAs and MaxAs note this function called using As 
threesigma = 3.*As .* DevAs;

if threesigma > 1
    threesigma = 1;     %If will cause negative size reduce integration range
end
MaxDp = 2.*(As + threesigma);
MinDp = 2.*(As - threesigma);

[ffa xxa AsCalc Density] = CalcSievedFormFunction(Model,MinDp,MaxDp,DevAs, f, c);

end

                
            
                