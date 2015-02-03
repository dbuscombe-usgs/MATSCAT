function [ffa xxa As Density] = CalcSievedFormFunction(Model,MinDp,MaxDp,DevDp, f, c)
%--------------------------------------------------------------------------
%   [ffa xxa Density] = CalcSievedFormFunction(Model,MinDp,MaxDp,DevDp, f,c)
% 
%   Calculates the required Form Function and scattering cross section for
%   a given frequency . The module takes into account the size distribution. 
%   The module will cope with an array of sediment size and frequency
%   Note that a normal size distribution is assumed
%--------------------------------------------------------------------------
%   Passed
%   Model = GLASS or SAND or USCGLASS 
%   MinDp = Small Sieve Diameter
%   MaxDp = Larger Sieve Diameter
%   DevDp = Deviation fraction to calculate standard deviation e.g. sigma=DevDp*MeanDp 
%   f = acoustic frequency in Hz [F1 F2 F3 ..]
%   c = speed of sound in water in m/s
%   Returns
%   ffa[NumSizes,NumFreq] = form function organised as ff[numSizes,numFreq]
%   xxa[NumSizes,NumFreq] = scattering cross-section as xx[numSizes,numFreq]
%   As = Mean Sediment radius (m)    
%   Density = Density of typical material
%-------------------------------------------------------------------------
%   Version
%   1.0     Initial Version
%   1.1     Modified for Toolkit to include lookup table for GLASS as was
%           very slow
%------------------------------------------------------------------------- 
    
    NumFractions = 50;  %Defines how many values to use in the calculation
    
%Ensure that the Frequency array is correct

[d1,d2] = size(f);
    if (d1 > 1)     %Matrix orientation wrong
       f = f';
       NumFreq = d1;
    else
       NumFreq = d2;
    end    

    %Ensure the MinDp and MaxDp array are the correct orientation
    [d1,d2] = size(MinDp);
    if (d2 > 1)     %Matrix orientation wrong
       MinDp = MinDp';
       MaxDp = MaxDp';
       NumSizes = d2;
    else
       NumSizes = d1;
    end    
 
    ffa = zeros(NumSizes,NumFreq);
    xxa = zeros(NumSizes,NumFreq);
    k = (2*pi*f/c);   %Has as many values as there are Freq
    
    if 1==strcmp('GLASS',Model) %determine min and max ka for lookup then create 101 lookup samples of ff xx
        for n = 1:NumSizes
            [f Dp Dpg] = Normal(MinDp(n),MaxDp(n),DevDp,NumFractions);
            ka_temp(n,:,:) = (Dp'/2)*k;
        end
        min_ka=min(min(min(ka_temp)));
        max_ka=max(max(max(ka_temp)));
        
        [ff_lookup xx_lookup] = TheoryFormFunction(min_ka:(max_ka-min_ka)/100:max_ka,c);
    end
    
    for n = 1:NumSizes
        [f Dp Dpg] = Normal(MinDp(n),MaxDp(n),DevDp,NumFractions);
        ka = (Dp'/2)*k;

        %Make suitable for multiple frequency calcs
        aa = repmat(Dp'/2,1,NumFreq);     
        f = repmat(f',1,NumFreq);
        %Now obtain the form function values for the given set of ka values
        if 1==strcmp('GLASS',Model)
            %[ff xx] = TheoryFormFunction(ka,c); %replaced from v1.0
            ff=interp1(min_ka:(max_ka-min_ka)/100:max_ka,ff_lookup,ka);         %interpolate lookups
            xx=interp1(min_ka:(max_ka-min_ka)/100:max_ka,xx_lookup,ka);
        elseif 1==strcmp('SAND',Model)
            [ff xx] = SedimentFormFunc(ka);
        end
        ffa(n,:) = (( sum(f.*aa) .* sum(f.*aa.^2 .* ff.^2) ) ./ sum(f.*aa.^3)).^0.5; 
        xxa(n,:) = (( sum(f.*aa) .* sum(aa.^2 .* xx) ) ./ sum(aa.^3)); 
        As(n) = Dpg/2;
    end
    %Example values returned user could user better estimates
     if 1==strcmp('GLASS',Model)
        Density = 2500;
     else
        Density = 2650; 
     end   
    
end

% Produces and array of values for the form funtion based on a given array
% of sediment radiuses
function [ff xx] = SedimentFormFunc(ka)
    x = ka;
    v1 = 0.25;
    x1 = 1.4;
    n1 = 0.5;
    v2 = 0.37;
    x2 = 2.8;
    n2 = 2.2;
    
    %Calculate Co
    Co = (1 - v1*exp( -((x - x1)/n1).^2)) .* ( 1 + v2*exp( -((x - x2)/n2).^2));
       
    Kf = 1.95; %From Haynes [1]
    Ka = 0.18; %From Haynes [1]
    ff = (Co * Kf .* x.^2) ./ (1 + Kf * x.^2);
    xx = ((4/3) * Ka * x.^4 ) ./ (1 + x.^2 + (4/3)*Ka.*x.^4);
    
end
    
function [f,Dp,Dpg]=Normal(r1,r2,DevDp,N)
% Input
% r2    MAX Diameter (in microns)
% r1    MIN Diameter (in microns)
% DevDp note sigma = DevDp*Dpg
% N     Number of fractions to use
% Output
% f     Frequency of occurrence of diameter Dp (size N)
% Dp    Diameter fraction with frequency f 
% Dpg   Mean value of population
% 
    Dpg = (r1 + r2)/2;  %Mean value of a Normal Distribution
 
    sg  = DevDp*Dpg;    %Standard deviation
  
    dp=[0:N-1];
    Dp = r1+(r2 - r1).*dp/N; %Assigned particle diameters
   
    %Equation for a normal distribution
    
    dNd = ( 1/(sg*(2*pi)^2)) .*exp(- ((Dp-Dpg).^2)/(2*sg.^2));
    
    f=(dNd)/sum(dNd); %Normalise

end


