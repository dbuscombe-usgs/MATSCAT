
function [SedAs SedMass Debug] = CalcSedimentSizeAndMass_efficient(Freq,Kt,r,V, StartBin, ff, XX, Density, As )
%*************************************************************************
%[SedSize SedMass] = CalcSedimentSizeAndMass(Model,DevAs,c,Freq,Kt,r,V,StartBin )
%
%Function to calculate the sediment size and mass from a range corrected
%profile containing multiple frequency data. 
%Note that this uses the equations for a 1/4phi size distributions
%Data should be preprocessed with RangeCorrectAbsProfiles(...)
%System information is held in the Abs structures such as
%the transducer frequencies, radius and Kt
%An important assumption is that all channels have the same set of bin
%ranges and therefore r is valid for all
%The dimenstions of the freq array must match that of the voltage array
%Passed:
%   Model       GLASS or SAND
%   DevAs       Estimated standard deviation for sample
%   c           Speed of Sound
%   Freq[Chan]
%   Kt[Chan]
%   r[Bin]      RangeArray (m) Assumes that it is speed of sound compensated 
%   V[Bin,Chan] Range Corrected Data organised as (Bin, Chan)
%   StartBin    The first bin to use (Best to start after 10cm)
%Returns   
%   SedAs[Bin]   Array of mean size (m)
%   SedMass[Bin] Array of mass concentration (g/l)
%   Debug        Structure containing the Debug information
% Version
%   1.0 Original Version
%   1.1 Added limits for As and Phi
%   1.2 Corrected problems with inversion equation and bin size integration
%   for attenuation
% ===============================================
% version 1.3
% modified Apr 2012 by Dan Buscombe so it calculates the form function externally - i.e.
% only once per burst rather than once per ping. It thereore now accepts
% the input arguments ff XX and Density
%=============================================
%*************************************************************************

% MinAs = 100E-6; %manually enter minimum and maximum mean radius size (in metres)
% MaxAs = 500E-6;

% PhiMin = -log (2000*MaxAs) / log (2);
% PhiMax = -log (2000*MinAs) / log (2);
% Phi = PhiMin:phi_increment:PhiMax;     %Generate Sediment sizes in fraction of Phi increments
% As = (2.^-Phi')/2000; %scale to radius in (m)
% 
% [ff XX Density] = CalcFormFunction(Model,As,DevAs,Freq,c);

NumChans = size(Freq,2);
NumBins  = size(r,1);
NumSizes = size(As,1);
BinSize = r(2) - r(1);
Ks = (ff./sqrt(repmat(As,1,NumChans)*Density));    %ff is 2D (As,Freq) => Make As 2D
AttenConstant = 3*XX./(4.*repmat(As,1,NumChans)*Density);   %Sediment Atten Constant

V = V ./repmat(Kt,NumBins,1);     % Apply Kt Correction to Voltages 
%Allocate Result Buffers
SedAs = zeros(NumBins,1);
SedMass = zeros(NumBins,1);
Debug.SedAtten = zeros(NumBins,NumChans);
Debug.SedDev = zeros(NumBins,NumSizes);

TotalAtten = zeros(1,NumChans);
atten = CalcBinSedimentSizeAndMass(StartBin,0); %Calc First Bin with no sediment attenuation correction
TotalAtten = r(StartBin)*atten;

for b = StartBin:NumBins
    atten = CalcBinSedimentSizeAndMass(b,TotalAtten);
    TotalAtten = TotalAtten + atten;
end

%Now use this value to recalculate the first bin having range corrected
%Then Repeat process for each bin passing the accumulated attenuation 

    function [atten] = CalcBinSedimentSizeAndMass(Bin,TotalAtten)
           %Calculate All Possible Mass Concentrations for Current Bin then
           %use to calculate attenuation to apply to next Bin
        M = (repmat( V(Bin,:) .*exp(4*TotalAtten),NumSizes ,1) ./ Ks).^2;
        
        s = std(M')./(mean(M'));    % Determine the Standard Deviation

        [~, I] = min(s); %min(s);            % determine the index of the most likely size distribution
        
        SedAs(Bin) = As(I);         %As radius

        SedMass(Bin) = mean (M(I,:));%Use Mean from all channels
        
        atten = .5*BinSize * AttenConstant(I,:) .* SedMass(Bin);

        %the following is provided for debuging to see how it processed
        Debug.As = As;                % Debug Array of Radius's tested
        Debug.SedAtten(Bin,:) = atten;  %Debug showing attenuation used
        Debug.SedDev(Bin,:) = s;    %Debug copy of deviation so can look at it
    end   
   
end
