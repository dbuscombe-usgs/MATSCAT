function [Vout r] = RangeCorrectAbsProfiles(T,S,D,r,Freq,At,Vin)
%*************************************************************************
% Function to apply a straight forward range correction that takes into 
% account the spreading loss and water attenuation. This helps identify the
% interesting sediment loadss
% Passed
%       T       Temperature degC
%       S       Salinity in pps
%       D       Depth in meters
%       r       Uncorrected Range
%       Freq    Channel Frequency
%       At      Effective Transducer Radius At
%       V       Array of Voltages from AQUAscat
% 
% Returns 
%       Vout    Range Corrected Voltage
%       r       Range Corrected for Speed of Sound
%*************************************************************************

    c = CalcSpeedOfSound(T,S,D);
    r = r * c/1500;     %Correct for the speed of sound
    
    % Calculate Water Attenuation
    atten = CalcWaterAttenuation(Freq,T,S,D);
        
    %Calculate the Nearfield correction
    Rn = (pi * (At.^2) .* Freq)/ c;
    z = r./Rn;
    Yeta = (1+1.35*z+((2.5*z).^3.2))./(1.35*z+((2.5*z).^3.2));
    Yeta = repmat(Yeta,1,size(Vin,2));
    r2 = repmat(r,1,size(Vin,2));
    Vout = Vin .* r2 .* Yeta .* exp(2*r2 * atten);
                    
end