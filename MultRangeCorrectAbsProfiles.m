function [V r] = MultRangeCorrectAbsProfiles(T,S,D,r,Freq,At,varargin)
%*************************************************************************
% Function to apply a straight forward range correction that takes into 
% account the spreading loss and water attenuation. This helps identify the
% interesting sediment loadss
% Passed
%       T       Temperature degC
%       S       Salinity in pps
%       D       Depth in meters
%       r       Uncorrected Range
%       Freq    Array of Frequencies
%       At      Array of Transducer Radius [At1 At2 etc]
%       Any number of 2D Voltage Arrays to average over
% 
% Returns 
%       Struct of V matrix's with corrected Voltage and Range
%*************************************************************************

        %Calculte the num of Profiles in the result
    c = CalcSpeedOfSound(T,S,D);
    r = r * c/1500;     %Correct for the speed of sound
    for n=1:size(varargin,2)
        V = varargin{n};
        atten = CalcWaterAttenuation(Freq(n),T,S,D);
      
        r2 = repmat(r,1,size(V,2));
        %Add the Nearfield correction
           % Deal with the Near Field Correction
        Rn = (pi * (At(n).^2) .* Freq(n))/ c;
        z = r./Rn;
        Yeta = (1+1.35*z+((2.5*z).^3.2))./(1.35*z+((2.5*z).^3.2));
        Yeta = repmat(Yeta,1,size(V,2));
        Vout = V .* r2 .* Yeta .* exp(2*r2 * atten);
        Ret(n).CorrectedV = Vout;
        Ret(n).CorrectedR = r;
    end
             
end