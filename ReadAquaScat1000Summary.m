function [ fOut ] = ReadAquaScat1000Summary(fIn , ExportType)
%--------------------------------------------------------------------------
% ReadAquaScat1000(fIn , ExportType)
% Import data from the AQUAscat1000 .aqa file format and writes it to
% a variety of Matlab format files. User specific file formats can be added
% If a base name is passed to the function then it is used otherwise the
% user is asked to select a file.
% The code is loosely based around code from CEV at UEA  
% The AQUAscat1000 file format start at 10 and packet formats at 4.
% A number of .MAT files are created to suit individuals, additional export
% methods could be applied. The type of export is controlled by the second
% parameter
%-------------------------------------------------------------------------
% Exports options:  
%   'SINGLE'   All Parameters and Data in a single <fIn>.mat file (Default)
%   'UEA'      Seperate File for parameters and individual frequencies
%   'STRUCT'   Generates the data using a structure
%   'USC'      Similar as in STRUCT but includes ABS.Info as used by the USC group
%-------------------------------------------------------------------------
% Notes
% ABS data from the AQUAscat organises its data slightly differently to the
% previous AQUAscat. It is now arranged in TimeSlots, thus if only one
% channel is enabled there will only be one timeslot. If two channels are
% enabled there will be two timeslots. The actual physical channel used
% during this time slot is recorded in the parameters. Note that it is
% possible to transmit on one channel and rx on another in some modes. 
%
% The checksum is now inverted compared to the original AQUAscat version
% and therefore not compatible if being tested. This change was made to
% improve searching for lost packets. Not implemented in this Import.
%-------------------------------------------------------------------------
% Uses:   
% ReadAquaScat1000(fIn,ExportType)-Define filename and Format
% ReadAquaScat1000(fIn)           -Define filename, use default Format
% ReadAquaScat1000(ExportType)    -Use GUI to get file, define Format
% ReadAquaScat1000                -USe GUI to get file, use default Format
%
Version = '1.4b';
%-------------------------------------------------------------------------
% VERSIONS
% 1.4b New file summary based on the ReadAquascat1000 function revision 1.4
%--------------------------------------------------------------------------
%
fprintf('\n-----------------------------------------------------------------');
fprintf('\nAQUAscat1000 File Summary Extraction Version %s by AQUATEC GROUP LIMITED',Version);
%
if nargin==0
    ExportType='SINGLE';
    [ftmp fpath] = uigetfile('*.aqa','Select Aquascat1000 Directory');
    
end
%   
if nargin==1
    if  (fIn(1:3)=='SIN' | fIn(1:3)=='UEA' | fIn(1:3)=='STR'| fIn(1:3)=='USC');
        ExportType=fIn;
        [ftmp fpath] = uigetfile('*.aqa','Select Aquascat1000 File');
      else
         ExportType='SINGLE';
         fname = [fIn,'.aqa'];      %assumes .AQA
         froot=fname;
    end
end
%    
if nargin==2
        fname = [fIn,'.aqa'];      %assumes .AQA
        froot=fname;
end
%
if nargin>2
    error(['Only Two Input Variables are allowed'])
end
%
fSearch = strcat(fpath,'*.aqa');
filesInDir= dir(fSearch);
filesSorted = sortrows({filesInDir.name}');
NumFiles = size(filesSorted,1);
fname = strcat(fpath,filesSorted{1});      %assumes .AQA
froot = strcat(fpath,regexprep (ftmp, '.aqa' ,'','ignorecase'));      %assumes .AQA


fid = fopen(fname,'rb');
if fid < 0 
    error(['Could Not Open File ',fname])
end
%
fprintf('\nFile being processed: %s',froot);
%
%Look for a File Version Packet, this is no longer right at the start as
%the logger personality is in that position. 
%
[status pktSize] = findAQUAscat1000Packet(fid,19);  %Read File Version from File 
if status == 1
    sdata=fread(fid,pktSize,'uint16');
    FileVersionMajor = sdata(2);
    FileVersionMinor = sdata(3);
    fprintf('\nAquaScat1000 File Version  %i.%i',sdata(2),sdata(3)); %GV-10
else
   fprintf(1,'%s','\nNo version packet therefore assume Version 5');
   FileVersionMajor = 5;
   FileVersionMinor = 0;  
end
clear sdata 
%
if 5 ~= FileVersionMajor
    fprintf('\n ***Not Version 5 - will probable fail***');
end
%
%--------------------------------------------------------------------------
%            Read in the Burst Start Time Information
%--------------------------------------------------------------------------
%
[status pktSize] = findAQUAscat1000Packet(fid,54);  %Read Burst Start Time 
if 1==status
    sdata=fread(fid,6,'uint16'); 
    WakeSource = fread(fid,1,'uint16');   %automatic wake up =3; ext.trigerred=6; GUI triger=4;
    BurstNumber = fread(fid,1,'uint32');
    Junk = fread(fid,1,'uint16');    
    BurstTime = datestr(sdata(1:6)');     %Time of first profile
    fprintf(1,'\nBurst Start Date: %s \nBurst Number: %i \nWake Source: %i \n',BurstTime,BurstNumber,WakeSource); %GV-10
else 
    BurstTime = datestr([0,0,0,0,0,0]);
    BurstNumber = 0;
    WakeSource = 0;
end
clear sdata;  
%
%--------------------------------------------------------------------------
%           Deal with reading the personality
%--------------------------------------------------------------------------
%
[status pktSize] = findAQUAscat1000Packet(fid,53);      %Read Personality Details
if 0==status
    fprintf('\nPersonality Packet not found abort');
    return
else
    PktStartPos = ftell(fid);
    skip = fread(fid,1,'uint16');   %Size of packet
    temp = fread(fid,1,'uint16');
    BoardVersion = sprintf('%i.%i',bitand(bitshift(temp,-8),255),bitand(temp,255));
    SerialNum =    sprintf('%i-%i',fread(fid,1,'uint16'),fread(fid,1,'uint16')); 
    MemorySize = fread(fid,1,'uint16'); 
    LoggerType = deblank(fread(fid,32,'uchar=>char'))';  %GV-10
    skip = fread(fid,1,'uint16');
    NumAbsChannels = fread(fid,1,'uint16');    % Number of ABS channels the system supports
    NumAuxChannels = fread(fid,1,'uint16');    % Number of AUX channels the system supports (8) 
    NumAbsTransducers = fread(fid,1,'uint16'); % Number of ABS Transducer Information that exist in personality table
    BatteryCapacity = fread(fid,1,'float');    %  not recorded at this time
    StandbyPower = fread(fid,1,'float');       %   "
    ActivePower = fread(fid,1,'float');        %   "
    %  None used space in here for future paramaters and licenses
    %  .....
    %  .....
    fseek(fid,PktStartPos + 112,'bof');     
    PtrToAuxInfo = fread(fid,1,'uint16');        % These are the offsets into the Packet
    PtrToTransducerInfo = fread(fid,1,'uint16'); %  
    % 
    % Read in the important information for the Aux Channels
    % First Need to assign the multiple dimension arrays
    %
    for i = 1 :NumAuxChannels 
        PtrToThisAux = PktStartPos + PtrToAuxInfo*2 + 400*(i-1);
        fseek(fid,PtrToThisAux,'bof');   %Move to the start of the ABS information
        skip = fread(fid,2,'uint16');
        AuxChannelName(i,:) = cleanString(fread(fid,16,'uchar =>char')); 
        skip = fread(fid,1,'uint16');
        AuxChannelUnit(i,:) = cleanString(fread(fid,8,'uchar =>char'))';
        skip = fread(fid,1,'uint16');
        AuxFlags(i) = fread(fid,1,'uint16');  
        skip = fread(fid,2,'uint16');
        AuxNumGain(i) = fread(fid,1,'uint16');
        skip = fread(fid,1,'uint16');
        AuxCalDate(:,i) = fread(fid,16,'char');
        skip = fread(fid,1,'uint16');
        AuxNumCoeff(i) = fread(fid,1,'uint16');
        skip = fread(fid,5,'uint16');
        
        fseek(fid,PtrToThisAux + 80,'bof'); %ensures aligned
        for j= 1 : AuxNumGain(i)            %Only import Coefficients for used gains
            AuxGainLabel(j,i,:) = cleanString(fread(fid,4,'uchar => char')); 
            skip = fread(fid,4,'uint16');
            AuxGainCoeff(:,j,i) = fread(fid,5,'float'); 
            AuxGainMin(j,i) =  fread(fid,1,'float'); % Minimum value used in calibration data
            AuxGainMax(j,i) =  fread(fid,1,'float'); % Maximum value used in calibration data            
            skip = fread(fid,10,'float'); 
        end
    end
    %
    %   Now Jump to the Transducer Info
    % 
    for i = 1: NumAbsTransducers
        fseek(fid,PktStartPos + PtrToTransducerInfo*2 + 200*(i-1),'bof');   %Move to the start of the ABS information
        TransducerSerialNum(i,:) = cleanString(fread(fid,20,'uchar => char')); % GV-10
        skip = fread(fid,1,'uint16');
        TransducerFrequency(i) = fread(fid,1,'float');  % In Hz
        TransducerRadius(i)    = fread(fid,1,'float');  % In meters
        TransducerBeamWidth(i) = fread(fid,1,'float');  % In Degrees (3dB beamdidth, derived from acoustic beam pattern)

        %
     %      1MHz   2MHz   2.5MHz  4MHz  5MHz
     % +/-  2.4    2.3    1.8     1.2   0.9  degrees beamwidth
     %      18mm   9.6mm  10mm    9.7mm  10.4mm Diameters 
     %
        skip = fread(fid,4,'float');
        TransducerKt(i) = fread(fid,1,'float');         % This is only if set in the personality
    end
      
end    
%
%--------------------------------------------------------------------------
%          Read in the Regime (Logger Set-Up) Information    
%--------------------------------------------------------------------------  
%
[status pktSize] = findAQUAscat1000Packet(fid,21);  %Read regime Details
if 0==status
    fprintf('\nRegime Packet not found abort');
    return
else
    PktStartPos = ftell(fid);
    %Session Information
    SessionControl = fread(fid,11,'uint16'); %Not interested in session start time at the moment%
    SessionTitle = cleanString(fread(fid,32,'uchar=>char'));
    AuxPingDiv = fread(fid,1,'uint16');      %The Aux channels are sampled at the PingRate divided by the auxPingDiv
    SerialPingDiv = fread(fid,1,'uint16');   %The serial pressure + Temperature sampling derived as above
    Junk = fread(fid,2,'uint16');
    NumAuxChans = fread(fid,1,'uint16');     %THIS IS INCORRECTLY SAVED BY THE SYSTEM, SHOULD BE 0 or 8
    NumAuxSamples = fread(fid,1,'uint32');
    
    % if NumAuxSamples~=0;     %  This is a temp trick to correct---------------------|
    %    NumAuxChans=8;       %  for NumAuxChans always being 0                      |
    % end                      %  It should be eliminated when corrected by Aquatec---|
    Junk = fread(fid,4,'uint16');
    PingRate = fread(fid,1,'uint16');        %This is the base profile rate before averaging
    NumPings = fread(fid,1,'uint32');        %This is number of profiles collected prior to averaging
    NumAbsTimeSlots = fread(fid,1,'uint16'); %Number of enabled ABS Channels
    Junk = fread(fid,3,'uint16');
    %
    PtrToAuxInfo = fread(fid,1,'uint16');    %These are the offsets into the Packet
    PtrToAbsInfo = fread(fid,1,'uint16');   
    %
    %  Calculate Aux Specific Information 
    %
    if 0 == AuxPingDiv          % if AuxPingDiv=0 then no aux channels are enabled and NumAuxChannels =0
        AuxSampleRate = 0;
        AuxNumSamples = 0;
        NumAuxChans=0;
    else 
        AuxSampleRate = PingRate / AuxPingDiv;
        AuxNumSamples = ceil(NumPings / AuxPingDiv); %GV-10
        if NumAuxChans ==0
            NumAuxChans=8;
        end
    end
    
    % Serial Pressure + Temperature information
    if 0==SerialPingDiv 
       NumSerialSamples = 0;
       SerialSampleRate = 0;        
    else
        NumSerialSamples = ceil(NumPings / SerialPingDiv);
        SerialSampleRate = PingRate / SerialPingDiv;
    end


    %Nothing useful in the channel
    %
    %   Now read in the ABS Channel information 
    %
    fseek(fid,PktStartPos + PtrToAbsInfo*2,'bof');      %Move to the start of the ABS information
   
    for j=1:NumAbsTimeSlots                             % For each ABS channel
        AbsComplex(j) = bitand(fread(fid,1,'uint16'),2);% For magnitude=0,complex=2
        AbsAverage(j) = fread(fid,1,'uint16');          % No of bursts averaged before saving
        AbsDecimation(j) = fread(fid,1,'uint16');       % Raw sampling rate along a profile is 19.2MHz, AbsDecimation i 
        AbsBinLengthMM(j) = 1.25*2^AbsDecimation(j);    % Converts to bin size in mm based on speed 1500ms-1
        AbsBinLength(j) = AbsBinLengthMM(j)/1500;       % Stored as time in seconds
        %
        %  Using the Trasnducer ID copy the relevent data across
        %
        TransducerId = fread(fid,1,'uint16')+1;     %Used to look up transducer information from personality
        AbsTransducerName(:,j)    = TransducerSerialNum(:,TransducerId);
        AbsTransducerRadius(j)    = TransducerRadius(TransducerId);    % in m
        AbsTransducerBeamWidth(j) = TransducerBeamWidth(TransducerId); % in degs
        AbsTransducerKt(j)        = TransducerKt(TransducerId);        
        %
        AbsTxFrequency(j)   = fread(fid,1,'float32');    %In Hz
        AbsRxFrequency(j)   = fread(fid,1,'float32');    %In Hz
        AbsTxPulseLength(j) = fread(fid,1,'float32');    %In seconds
        Junk = fread(fid,1,'float32');
        AbsStartingGain(j)  = fread(fid,1,'float32');    %In dB with reference to default (built-in) Gain of system
        AbsTVG(j) = fread(fid,1,'float');                %In dB / bin where first bin has StartingGain (not used, =0)
        powerlevel=fread(fid,1,'uint16');
        AbsPowerLevelpc(j) = 100/2^(2*powerlevel);       %Power Level in % of Maximum
        AbsPowerLevel(j) = -20*log10(2^powerlevel);      %Power Level in dB relative to Maximum Power
        AbsStartBin(j)   = fread(fid,1,'uint16');        %Number of Bins from Start of Tx pulse before recording
        AbsNumBins(j)    = fread(fid,1,'uint16');        %Number of Bins recorded
        AbsRxChan(j)     = fread(fid,1,'uint8');         %Indicates which channel on AQUAscat used for TX
        AbsTxChan(j)     = fread(fid,1,'uint8');         %Indicates which channel on AQUAscat used for RX
        %
        Junk = fread(fid,12,'uint16');
        %
        AbsNumProfiles(j) = NumPings / AbsAverage(j);    %Calculate the number of profiles that should be recorded for this channel
        AbsProfileRate(j) = PingRate / AbsAverage(j);    %Calculate the stored profile rate           
              
  %      if AbsTxPulseLength(j) ==0 
  %              AbsTxPulseLength(j) = 2*AbsBinLengthMM(j)/(1000*1500);  %t =(2-way travel distance)/c
  %      end
         AbsBinRange(:,j) = [AbsStartBin(j):1:AbsStartBin(j)+AbsNumBins(j)-1]' *AbsBinLengthMM(j)/1000; %in m --

    end   
end


%--------------------------------------------------------------------------
% Now open each *.aqa file in the directory and summarise the data
% This function assumes full burst length
%--------------------------------------------------------------------------

fclose(fid);
AbsMean = zeros(AbsNumBins(1),NumFiles,NumAbsTimeSlots);
AbsMin = zeros(AbsNumBins(1),NumFiles,NumAbsTimeSlots);
AbsMax = zeros(AbsNumBins(1),NumFiles,NumAbsTimeSlots);
AbsStd = zeros(AbsNumBins(1),NumFiles,NumAbsTimeSlots);
AuxMean = zeros(NumFiles,NumAuxChans);
PressMean = zeros(NumFiles,2);

for n = 1 : NumFiles
    [AbsTmp AuxTmp PressTmp] = ReadAQUAscat1000DataOnly(filesSorted{n});
    for i = 1:NumAbsTimeSlots
        AbsMean(:,n,i) = mean(AbsTmp(:,:,i)')'; 
        AbsMin(:,n,i) = min(AbsTmp(:,:,i)')'; 
        AbsMax(:,n,i) = max(AbsTmp(:,:,i)')'; 
        AbsStd(:,n,i) = std(AbsTmp(:,:,i)')'; 
    end
if NumAuxSamples > 0 
    AuxMean(n,:) = mean(AuxTmp);
end
if NumSerialSamples > 0
    PressMean(n,:) = mean(PressTmp);
end
end

froot = strcat(fpath,'Summary');
if strcmp(ExportType,'SINGLE')
        fOut = ExportSingle();
elseif strcmp(ExportType,'UEA')
        fOut = ExportUEA();
elseif strcmp(ExportType,'STRUCT')
        fOut = ExportStruct();
elseif strcmp(ExportType,'USC')  
        fOut = ExportUSC();
end

fprintf('\nFinished AQUAscat Import');
fprintf('\n-----------------------------------------------------------\n');

       



%--------------------------------------------------------------------------
%   Function to deal with reading in all the data from a specific File
%   This assumes knowledge of the Data content from a previous read.  
%   Returns the entire data block
%--------------------------------------------------------------------------

function [AbsData AuxData PressTempData ] = ReadAQUAscat1000DataOnly(ftmp)
    fname = strcat(fpath,ftmp);      %assumes .AQA
    froot = strcat(fpath,regexprep (ftmp, '.aqa' ,'','ignorecase'));      %assumes .AQA
    fid = fopen(fname,'rb');
    if fid < 0 
        error(['Could Not Open File ',fname])
    end
    %
    fprintf('\nFile being processed: %s',froot);
    

    %  Allocate Memory for the Data
    AuxData = zeros(AuxNumSamples,NumAuxChans);
    AbsData = zeros(AbsNumBins(1),AbsNumProfiles(1),NumAbsTimeSlots);
    PressTempData = zeros(NumSerialSamples,2);

    AuxIndex = 0;
    SerIndex = 0;
    AbsIndex = zeros(NumAbsTimeSlots,1);
    fseek(fid,0,'bof');
    %
    %   Now Read in all the Data
    %  
    while(0 == feof(fid)) 
        [status pktType pktSize] = readNextAQUAscat1000Header(fid);
        if 1 == status
            switch(pktType)
                %
                % Case 22: Data were saved as magnitude values (normalize by 65536)
                %
                case (22)
                    chan = fread(fid,1,'uint16') + 1;
                    AbsIndex(chan) = AbsIndex(chan) + 1;    %Increase the Index
                    AbsData(:,AbsIndex(chan),chan) = fread(fid,AbsNumBins(chan),'uint16')/65536;
                %
                % Case 41: Data were saved as Complex values (normalize by 32768)
                %
                case (41)
                    chan = fread(fid,1,'uint16') + 1;
                    AbsIndex(chan) = AbsIndex(chan) + 1;    %Increase the Index
                    sss = fread(fid,2*AbsNumBins(chan),'int16')/32768;
                    AbsData(:,AbsIndex(chan),chan) = sss(1:2:end)+sqrt(-1)*sss(2:2:end);
                %
                % Case 46: Auxiliary Channel Data
                %
                case (46)
                    temp = fread(fid,1,'uint16');   %Gain settings
                    AuxGain = [bitand(temp,3)+1,bitand(bitshift(temp,-2),3)+1,bitand(bitshift(temp,-4),3)+1,bitand(bitshift(temp,-6),3)+1,1,1,1,1];
                    Junk = fread(fid,1,'uint16');   %Flags
                    AuxIndex = AuxIndex + 1;
                    AuxData(AuxIndex,1:NumAuxChans) = fread(fid,NumAuxChans,'uint16')'; 
                %   
                %
                % Case 55: External Pressure Channel for upgraded ABS system -
                % details to be provided
                %
                case (55)
                    chan = fread(fid,1,'uint16') + 1;
                    SerIndex = SerIndex + 1;    %Increase the Index
                    PressTempData(SerIndex,:) = fread(fid,2,'float')';

                otherwise
                    fseek(fid,pktSize*2,'cof');
            end
        end
    end
    %
     %
    % Need to apply calibration coefficients for the Aux Channels, now that the
    % gain is fixed, auto gain is not supported. Therefore will use last value
    %
    for i=1:NumAuxChans
        coeff = AuxGainCoeff(:,AuxGain(i),i);  %Matlab orders polynomial coeff opposite to Aquatec
        coeff = reshape(coeff,1,5);
        coeff = fliplr(coeff);
        AuxData(:,i) = polyval(coeff,AuxData(:,i));
    end
    fprintf(' Num Abs: %i  Num Aux: %i  Num Serial %i',AbsIndex(1),AuxIndex,SerIndex);    
    fclose(fid);    % Close the *.aqa data file
end




%--------------------------------------------------------------------------
%            NOW DEAL WITH THE EXPORTING IN THE DIFFERENT FORMATS
%---------------------------------------------------------------------------
%
    function [fName]= ExportSingle()
        fprintf('\nSaving Data using SINGLE Filter');
        %Deal with Saving all the information from this file
        fName =strcat(froot,'.mat'); %save everything in a single MAT file
        eval(['save ' fName ' SessionTitle  BurstTime BurstNumber WakeSource PingRate NumPings '...
         'NumAuxChans NumAuxSamples AuxSampleRate AuxNumSamples '...
         'AuxChannelName AuxChannelUnit AuxData '...
         'NumAbsTimeSlots AbsNumProfiles AbsComplex AbsAverage AbsProfileRate '...
         'AbsDecimation AbsBinLengthMM AbsBinLength AbsStartBin AbsNumBins '...
         'AbsTransducerName AbsTransducerRadius AbsTransducerBeamWidth AbsTransducerKt ' ...
         'AbsTxFrequency AbsRxFrequency AbsTxPulseLength '...
         'AbsStartingGain AbsTVG AbsPowerLevel  AbsRxChan AbsTxChan '...
         'AbsMean AbsData AbsBinRange']);
     
    end
%--------------------------------------------------------------------------
% THE EXPORT BELOW (UEA)IS DESIGNED TO PROVIDE THE DATA IN THE SAME FORMAT AS 
% USED BY UEA FOR THEIR PROCESSING OF OLD AQUASCAT FILES. NOTE THAT CERTAIN
% VARIABLES NO LONGER RELEVANT. THIS EXPORTS THE VARIABLES INTO SEPERATE 
% PARAMETER AND FREQUENCY FILES
%--------------------------------------------------------------------------
function [fmeanfile] = ExportUEA()
    fprintf('\nSaving Data using UEA Filter');
   PingRate = PingRate;
   Averaging = AbsAverage(1);           %Assume all channels have same averaging 
   ControlFlag1 = 0;
   ControlFlag2 = 0;
   ProfileRecordRate = PingRate / Averaging;
   NumOfChannels = NumAbsTimeSlots;
   TXPulseLength = AbsTxPulseLength;    % NOW IN SECONDS DIVIDE BY 
   Decimation = AbsDecimation;
   BinSize = AbsBinLengthMM;
   StartRange = AbsStartBin .* AbsBinLengthMM / 10; % Start distance in CM 
   StopRange = (AbsStartBin + AbsNumBins).*AbsBinLengthMM / 10;
   BinsPerProfile = AbsNumBins;
   PulseInterval = 1 ./AbsProfileRate;
   PulseFrequency = AbsProfileRate;
   TransAngle = 0;
   ScaledGain = 0;  %Not Relevent
   GainStepPerCm = 0; %Not relevent
   FixedGain = AbsStartingGain;
   GainStep = AbsTVG;   %Gain change per bin
   NumberOfProfiles = AbsNumProfiles;
   
   [yyear mmonth dday hhour mminute ssecond] = datevec(BurstTime);
   
   fparam =strcat(froot,'_Param.mat'); %save the parameters to go with the data
   save (fparam, 'yyear', 'mmonth','dday','hhour','mminute','ssecond','PingRate',...
        'Averaging','ControlFlag1','ControlFlag2','ProfileRecordRate',...
        'NumOfChannels','TXPulseLength','Decimation','BinSize','StartRange','StopRange',...
        'BinsPerProfile','PulseInterval','PulseFrequency','TransAngle','ScaledGain',...
        'GainStepPerCm','FixedGain','GainStep','NumberOfProfiles')
  %Option To allow Data to be output in seperate files 
  if 1 <= NumAbsTimeSlots
      F_1 = AbsMean(:,:,1);
      absFilename = sprintf('%s_F1.mat',froot);
      save(absFilename, 'F_1');
  end
  if 2 <= NumAbsTimeSlots
      F_2 = AbsMean(:,:,2);
      absFilename = sprintf('%s_F2.mat',froot);
      save (absFilename, 'F_2');
  end
  if 3 <= NumAbsTimeSlots
      F_3 = AbsMean(:,:,3);
      absFilename = sprintf('%s_F3.mat',froot);
      save (absFilename, 'F_3');
  end
  if 4 <= NumAbsTimeSlots
      F_4 = AbsMean(:,:,4);
      absFilename = sprintf('%s_F4.mat',froot);
      save ( absFilename, 'F_4');
  end
   
    fmeanfile=strcat(froot,'_fmean');
    save ( fmeanfile, 'AbsMean');

end

%--------------------------------------------------------------------------
% This Export Format (STRUCT) saves the info in structure format
%--------------------------------------------------------------------------
function [fName] = ExportStruct()
    fprintf('\nSaving Data using STRUCT Filter');
    %Deal with Saving all the information from this file
    fName =strcat(froot,'.mat'); %save everything in a single MAT file
  
    %Generate a structure for each ABS channel
    for i = 1:NumAbsTimeSlots
       %Abs(i).NumProfiles = AbsNumProfiles(i); 
       %Abs(i).ProfileRate = AbsProfileRate(i);
       Abs(i).NumBursts = NumFiles;
       Abs(i).Freq = AbsTxFrequency(i);
       Abs(i).At = AbsTransducerRadius(i);
       Abs(i).MeanV =AbsMean(:,:,i);
       Abs(i).MinV =AbsMin(:,:,i);
       Abs(i).MaxV =AbsMax(:,:,i);
       Abs(i).StdV =AbsStd(:,:,i);
       Abs(i).r = [AbsStartBin(i):1:AbsStartBin(i)+AbsNumBins(i)-1]' *AbsBinLengthMM(i)/1000; %in m
       Abs(i).TxPulseLength = AbsTxPulseLength(i);
       Abs(i).RxGain = AbsStartingGain(i);    %in dB
       Abs(i).TxGain = AbsPowerLevel(i);  % as 
    end
    %Generate a structure for each Auxiliary channel
    if NumAuxChans==8
        Aux.Temperature = AuxData(:,5);
        Aux.Battery     = AuxData(:,8);
        Aux.Pressure    = AuxData(:,6);
        Aux.Analogue01  = AuxData(:,1);
        Aux.Analogue02  = AuxData(:,2);
        Aux.SampleRate  = AuxSampleRate;
        save(fName, 'Abs','Aux');
    else
        save(fName, 'Abs');
    end
end
%--------------------------------------------------------------------------
% Other Export Formats
%--------------------------------------------------------------------------
function [fName] = ExportUSC()
    fprintf('\nSaving Data using USC Filter');
    %Deal with Saving all the information from this file
    fName =strcat(froot,'.mat'); %save everything in a single MAT file
  
    %Generate a structure for each ABS channel
    for i = 1:NumAbsTimeSlots
       Abs(i).NumBursts = NumFiles; 
  %     Abs(i).ProfileRate = AbsProfileRate(i);
       Abs(i).Freq = AbsTxFrequency(i); %---
       Abs(i).At = AbsTransducerRadius(i); % --
       Abs(i).MeanV =AbsMean(:,:,i);
       Abs(i).MinV =AbsMin(:,:,i);
       Abs(i).MaxV =AbsMax(:,:,i);
       Abs(i).StdV =AbsStd(:,:,i);
       Abs(i).r = [AbsStartBin(i):1:AbsStartBin(i)+AbsNumBins(i)-1]' *AbsBinLengthMM(i)/1000; %in m --
       Abs(i).TxPulseLength = AbsTxPulseLength(i);
       Abs(i).RxGain = AbsStartingGain(i);    %in dB  ---
       Abs(i).TxGain = AbsPowerLevel(i);  % in dB of Transmit (0 =100% power) 
    end
    %
    % Generate the ABS_Info structure
    %
    ABS_Info.BinSize_mm=AbsBinLengthMM;
    ABS_Info.Fixed_Gain_db=AbsStartingGain;
    ABS_Info.Gain_Step_db=AbsTVG;   %In dB / bin where first bin has StartingGain
    ABS_Info.AbsPowerLevel= AbsPowerLevel;  %Power Level in % of Maximum
    ABS_Info.No_of_Freqs=NumAbsTimeSlots;
    ABS_Info.StartRange_cm=AbsStartBin.*AbsBinLengthMM/10;
    ABS_Info.StopRange_cm=(AbsStartBin(i)+AbsNumBins(i)-1)*AbsBinLengthMM/10;
    ABS_Info.Transducer_diameter_mm=AbsTransducerRadius*1000;  % convert from m to mm
    ABS_Info.Transmit_Frequency_MHz=AbsTxFrequency/1000000; % convert from Hz to MHz
    %
    %Generate a structure for each Auxiliary channel
    %
    % I MIGHT HAVE TO CHANGE THIS FOR THE UPGRADED SYSTEM AS PRESSURE MIGHT
    % ADD ANOTHER CHANNEL
    %
    if NumAuxChans==8
        Aux.Temperature = AuxData(:,5);
        Aux.Battery     = AuxData(:,8);
        Aux.Pressure    = AuxData(:,6);
        Aux.Analogue01  = AuxData(:,1);
        Aux.Analogue02  = AuxData(:,2);
        Aux.SampleRate  = AuxSampleRate;
        save(fName, 'Abs','Aux','ABS_Info');
    else
        save(fName, 'Abs','ABS_Info');
    end
    

end

    
%
%--------------------------------------------------------------------------
% The following functions are used by the main function to search and read
% file packet headers. 
%--------------------------------------------------------------------------
function [Status PktType PktSize] = readNextAQUAscat1000Header(fid)   
PktType=fread(fid,1,'uint8');      %
PktVersion=fread(fid,1,'uint8');    %Packet Version
PktSize=fread(fid,1,'uint16'); %size of packet body
PktChecksum=fread(fid,1,'uint16');
if 0 == feof(fid)
    Status = 1;
else
    Status = 0;
    %Need to test the checksum
end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [status pktSize] = findAQUAscat1000Packet(fid, type)
%
% This function looks for the specific packet type and returns once it has
% been found. It always starts from the beginning of the file and reads in
% a header at a time. If header not required then jumps to next.
%
fseek(fid,0,'bof');
while (0 == feof(fid))
    [status pktType pktSize] = readNextAQUAscat1000Header(fid);
    if  0 == status
        break;
    else if pktType == type 
            break;
        else
            fseek(fid,pktSize*2,'cof'); %Skip past packet content
        end
    end
end
end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function [strOut] = cleanString(strIn)
strOut = char(zeros(size(strIn))');
strTmp = deblank(strIn);
strOut(1:size(strTmp)) = strTmp;
end 
end
  