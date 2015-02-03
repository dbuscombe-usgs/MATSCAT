function [ff xx] = CalcSedimentFormFunction(As, f, c)
%--------------------------------------------------------------------------
%   Calculates the required Form Function and scattering cross section for
%   a given frequency . The module takes
%   into account the size distribution. 
%   The module will cope with an array of sediment size
%--------------------------------------------------------------------------
%   Passed
%   As = Mean partical Radius (m) (As1 As2 As3)
%   f = acoustic frequency in Hz [F1 F2 F3 ..]
%   c = speed of sound in water in m/s
%   Returns
%   ff = form function organised as ff[numSizes,numFreq]
%   xx = scattering cross-section as xx[numSizes,numFreq]
%-------------------------------------------------------------------------
%   Version
%   1.0     Initial Version
%-------------------------------------------------------------------------
%   References
%  [1] E.D. Thoeston and D.M. Hanes, J.A.S.A, 104(2) Pt1 1998
%

    NumFreq = size(f,2);
    NumSizes = size(As,1);
    
    k = 2*pi*f/c;   %Has as many values as there are Freq
    x = repmat(k,NumSizes,1) .* repmat(As,1,NumFreq);     

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


