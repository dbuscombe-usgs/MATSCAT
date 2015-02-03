
function [] = range_correct_abs(basepath,temp,depth,salinity,AbsKt)

% range_correct_abs.m
% A function to take raw acoustic amplitudes and return 
% range-corrected acoustic backscatter amplitudes,
%
% INPUTS:
% 1. basepath = a string containing the directory of files which you would
% like to process (or a cellular array of strings with the names of
% multiple directories - see batch_range_correct_abs.m)
% 2. temp = scalar or cellular array of scalars. water temperature in degrees C
% 3. depth = scalar or cellular array of scalars. water depth in metres
% 4. salinity = scalar or cellular array of scalars. in pss
% 5. AbsKt = a numeric array (or cellular array of numeric arrays if used in batch mode) 
% of acoustic attenuation coefficient per transducer, as provided by the instrument manufacturer
%
% OUTPUTS:
% A structure called 'data' is written to a file the same as the input
% file appended by "_rangecorrectABS.mat" in a subdirectory 'data_out/processed/':
% the structure contains the following fields:
%                 T: temperature in degrees C (same as input)
%                 D: water depth in metres (same as input)
%                 S: salinity in PSS (same as input)
%             AbsKt: 1xnumber of ttransducers array of transducer coefficients (same as input)
%     start_mattime: 1x1 time of first ping in matlab datenum format
%      start_string: string of time of first ping
%           mattime: timestamp per ping in matlab datenum format
%               Abs: structure containing instrument info:
%                 a 1xnumber of transducers struct array with fields:
%                     NumProfiles
%                     ProfileRate
%                     Freq
%                     At
%                     V
%                     MeanV
%                     r
%                     TxPulseLength
%                     RxGain
%                     TxGain
%                     Average
%             Vrcor: [number of bins x number of pings x number of
%             transducers] array of acoustic backscatters
%             rrcor: [number of bins x number of transducers] array of
%             distances
%                 c: calculated speed of sound
%              ntim: 1x1 number of pings 
%
% Daniel Buscombe May-June 2012

for f=1:length(basepath)

    files=ReadImDir(basepath{f},'aqa');
    addpath(basepath{f})

    mkdir([basepath{f},filesep,'data_out',filesep,'raw'])
    mkdir([basepath{f},filesep,'data_out',filesep,'processed'])

    for i=1:size(files,1)

        data.T=temp{f};
        data.D=depth{f};
        data.S=salinity{f};

        data.AbsKt(1) = AbsKt{f}(1);
        data.AbsKt(2) = AbsKt{f}(2);
        data.AbsKt(3) = AbsKt{f}(3);

        try
        ReadAquaScat1000([files(i,1:end-4)],'STRUCT',[basepath{f},filesep,'data_out',filesep,'raw'])

        load(deblank([basepath{f},filesep,'data_out',filesep,'raw',filesep,files(i,:),'.mat'])) 

        tstring=files(i,1:end-4);
        [data.start_mattime,data.start_string]=get_abs_timestamp(tstring);

        time=datenum(BurstInfo.StartTime):1/24/60/60/Abs(1).ProfileRate:datenum(BurstInfo.StartTime)+1/24;

        data.mattime=interp1([1:length(time)],time,[1:length(Abs(1).V)]);

        data.Abs=Abs;

        index_use=find(mean(Abs(1).V(:,1:end))>0);

        ndim=size(Abs,2);%This is probably always three but...

        data.Vrcor=zeros(length(Abs(1).r),length(index_use),ndim);
        data.rrcor=zeros(length(Abs(1).r),ndim);
        for k=1:ndim
            [data.Vrcor(:,:,k),data.rrcor(:,k)] = RangeCorrectAbsProfiles(data.T,data.S,data.D,data.Abs(k).r,...
                data.Abs(k).Freq,data.Abs(k).At,data.Abs(k).V(:,index_use));
        end

        data.c = CalcSpeedOfSound(data.T,data.S,data.D);
        data.ntim=size(data.Vrcor,2);

        save([basepath{f},filesep,'data_out',filesep,'processed',filesep,files(i,:),'_rangecorrectABS.mat'],'data')

        catch
            continue
        end

    end

end





