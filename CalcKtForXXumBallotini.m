% Calibration Script for calculating Kt using a near homogenous 1/4phi size
% distributed ballotini. 
% Note that these functions use the AQUAscat data extracted in the 'STRUCT' format
% It is indended that the user will copy the script and alter the fields
% marked ***field*** to match their calibration. 

clear
close all
% Load the background data and calculate the Mean
load('DemoBackground.mat')          % ***Enter Actual Background Filename For Test***
BackgroundAbs = Abs;

% Load the Calibration File
load('DemoCalFile.mat')      % ***Enter Actual Calibration Filename For Test***



NumChans = size(Abs,2)
NumBins = size(Abs(1).r,1)

Frequencies = [Abs.Freq]


% Sediment Details including reference Range
RefRange = 0.7;                 % *** Range at which the pump samples were taken (m)***
RefMass = 0.281;                % *** Result of pump samples (g/l) *** 
RefAverageBins = 20;            % *** Number of bins to average over either side of Reference ***
As = 137E-6                     % *** Mean Sediment Radius (m)
MinAs = As *0.91;               % Min and Max size of uniform distribution  (assume 1/4 phi)
MaxAs = As *1.09;

% Water Details
T = 15.5                        % *** Water Temperature degC ***
S = 0;                           % Assume Zero Salinity (ppt)
D = 0;                           % Assume Zero Depth (m) 

% Check the voltage levels and provide the user with an indication 
% of the range that channels should no longer be used.
% The bins effected are stored in a variable called InvalidBins so that
% they can be used at a later stage to mask the Kt values

WarningLevel = 20/65535;        % If voltage falls below 20 ADC counts warn user 
for ch = 1:4
    index = find(Abs(ch).MeanV < WarningLevel);
    if (~isempty(index))
       fprintf('\nChannel %i signal too small after %0.2f m',ch,Abs(ch).r(index(1))); 
       Abs(ch).ValidToBin = index(1);
       Abs(ch).InvalidBins = index;
       %Strip out data that is below the warning level
    else
        Abs(ch).ValidToBin = size(Abs(ch).r,2);
    end
    
    index = [];
end
    
% The following removes the background noise from the data and then
% corrects the data for the near field and spherical spreading as well as
% Water attenuation and the spead of sound
% Note that only the index's where the array is below the warning threshold are used. 
for ch = 1:4
    V = real((Abs(ch).MeanV .^2 - BackgroundAbs(ch).MeanV .^2).^0.5); %Take off background
    [Abs(ch).CorrectedV Abs(ch).r] = RangeCorrectAbsProfiles(T,S,D,Abs(ch).r,Abs(ch).Freq,Abs(ch).At,V);
end

%---- Plot Range Corrected Voltage ----
h = figure;
subplot(1,4,1)
plot([Abs.CorrectedV],[Abs.r])
set(gca,'YDir','reverse'); 
v = axis;
axis([0 0.01 0 v(4)]) ;
xlabel('Range Corrected Voltage')
ylabel('Range (m)')
%--------------------------------------


%Calculate Useful Information
   BinSize = Abs(1).r(2) - Abs(1).r(1);
   %Determine which bins to extract the Kt value from
   RefBin = int16((RefRange - Abs(ch).r(1)) / (Abs(ch).r(2) - Abs(ch).r(1)));
   RefIndex  = RefBin-RefAverageBins:RefBin+RefAverageBins;

for ch=1:4
    [ff XX Density] = CalcBallotiniFormFunction(MinAs,MaxAs,Abs(ch).Freq, CalcSpeedOfSound(T,S,D));
    Abs(ch).SedConstant = (3*XX)/(4*Density*As); %In Neper/m
    Abs(ch).Ks = ff/sqrt(As*Density);
    
end

%Now calculate Kt assuming constant mass concentration
for ch=1:4
    AttenSed = RefMass*Abs(ch).SedConstant; %In Neper/m
    Abs(ch).ArrayKt = (Abs(ch).CorrectedV .* exp(2*Abs(ch).r*AttenSed))./(Abs(ch).Ks*sqrt(RefMass));
    Abs(ch).ArrayKt(Abs(ch).InvalidBins) = 0;
    Abs(ch).Kt = mean(Abs(ch).ArrayKt(RefIndex));     %average ±5bins from reference range
end 

%----Plot Kt-----
subplot(1,4,2)
plot([Abs.ArrayKt],[Abs.r])
set(gca,'YDir','reverse'); 
v = axis;
axis([0 0.02 0 v(4)]);
xlabel('System Constant Kt (Uncorrected)')
ylabel('Range (m)')
%------------------

%Now Calculate Mass using Kt and one interation
SedConstants = repmat([Abs.SedConstant],NumBins,1);
M = ( [Abs.CorrectedV] ./ repmat([Abs.Ks].*[Abs.Kt],NumBins,1)).^2; %First assume SedAtten =0
M = M .* exp (4 * BinSize * cumtrapz( SedConstants .*M));

%----Plot Mass-----
subplot(1,4,3)
plot(M(:,1),[Abs(ch).r])
%plot(M,[Abs.r])        %Use this to plot all Mass Results
set(gca,'YDir','reverse'); 
v = axis;
axis([0 1 0 v(4)]);
xlabel('Mass Concentration g/l')
ylabel('Range (m)')
%------------------

%Now rework out Kt for all channels using Mass from Channel 1
M = repmat(M(:,1),1,NumChans);
Atten = BinSize * cumtrapz( SedConstants.*M);

for ch=1:4
    Abs(ch).ArrayKt = (Abs(ch).CorrectedV .* exp(2 * Atten(:,ch)))./(Abs(ch).Ks*(M(:,ch).^0.5));
    Abs(ch).ArrayKt(Abs(ch).InvalidBins) = 0;
    Abs(ch).Kt = mean(Abs(ch).ArrayKt(RefIndex));     %average ±5bins from reference range
end
 
KT = [Abs.Kt]

%----Plot Kt-----
subplot(1,4,4)
plot([Abs.ArrayKt],[Abs.r])
set(gca,'YDir','reverse'); 
v = axis;
axis([0 0.02 0 v(4)]);
xlabel('System Constant Kt (Corrected with 1MHz)')
ylabel('Range (m)')
%------------------

%*************************************************************************
% Now create a new graph to demonstrate calculating the Mass using the Kt 
%*************************************************************************
figure


%Ignore the first 20 bins as sediment variable
%BinsToUse = 10:size(Abs(1).r,1);
%for ch = 1:4
%    Abs(ch).r = Abs(ch).r(BinsToUse);
%    Abs(ch).CorrectedV =     Abs(ch).CorrectedV(BinsToUse);
%end
%NumBins = size(Abs(1).r,1)

%Calculate Mass first assuming assume no sediment attenuation
M = ( [Abs.CorrectedV] ./ repmat([Abs.Ks].*[Abs.Kt],NumBins,1)).^2;

%---- Plot First ----
subplot(1,3,1);
plot(M,[Abs.r]);
set(gca,'YDir','reverse'); 
xlabel('Concentration (Iteration)')
ylabel('Range (m)')
hold on;

subplot(1,3,2);
%plot(M,[Abs.r]);
set(gca,'YDir','reverse'); 
xlabel('Sediment Attenuation')
ylabel('Range (m)')
hold on;
%--------------------

LastAtten = zeros(NumBins,NumChans);
r = [Abs.r];
SedConstants = repmat([Abs.SedConstant],NumBins,1);
%Now using the result for M calculate the attenuation and repeat
for i=1:10
    Atten = BinSize * cumtrapz(SedConstants .*M);    
    %Atten(1,:) = Atten(1,:) .* r(1,:);   %Use value in first bin to correct blanking distance
    %To avoid problems with first 20 bins corrupting the data
   %Atten(1:40,:) = repmat(Atten(40,:),40,1);
    M = M .* exp( 4*(Atten - LastAtten));
    LastAtten = Atten;
    subplot(1,3,1);
    plot(M,[Abs.r]);
    subplot(1,3,2);
    plot(Atten,[Abs.r]);

end

% Set The Axis Scaling
subplot(1,3,1);
v = axis;
axis([0 0.5 0 v(4)]);

% Set The Axis Scaling
subplot(1,3,2);
v = axis;
axis([0 1 0 v(4)]);

%---Now Plot the final itteration ---
subplot(1,3,3)
plot(M,[Abs.r])
set(gca,'YDir','reverse'); 
xlabel('Mass Concentration (g/l)')
ylabel('Range (m)')
v = axis;
axis([0 0.5 0 v(4)]);





