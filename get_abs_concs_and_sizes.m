
function []=get_abs_concs_and_sizes(basepath,usefiles,Params,StartBin,phi_increment,Model)

% get_abs_concs_and_sizes.m
% A function to take range-corrected acoustic backscatter amplitudes,
% calculate a form function based on user inputs, and return calculated
% profiles of sediment mass concentration and mean size
%
% INPUTS:
% 1. basepath = a string containing the directory of files which you would
% like to process (or a cellular array of strings with the names of
% multiple directories - see batch_get_abs_concs_and_sizes.m)
% 2. usefiles = a numeric array of indexs of files in the directory specified
% by basepath (or a cell of numeric indices if basepath is a cellular array
% - see batch_get_abs_concs_and_sizes.m) which you wish to process
% 3. Params = a numeric array (or cellular array of numeric arrays if used
% in batch mode - see batch_get_abs_concs_and_sizes.m) which contains a vector of 3 numbers: 
% StDev, min grain size, and max grain size (in metres)  
% 4. StartBin = scalar indicating which bin (increasing from the
% transducer) to start calculations (usually the top few bins are
% unreliable)
% 5. phi_increment = increment (in Phi units (-log base 2 of mm)) between
% min and max grain sizes at which the sizes are calculated
% 6. Model = either 'GLASS' or 'SAND' depending on the substrate. Note
% these functions have not been tested for silts/muds
%
% OUTPUTS:
% The following variables are written to a file the same as the input
% file appended by "_size_n_mass.mat":
% 'Size','Mass','x','StartBin','Params','Model','phi_increment'
% Size = particle diameter in mm (mxn, m=bins from transducer (bottom to top), t=pings
% (time))
% Mass = mass concentration (grams per litre) (mxn, m=bins from transducer (bottom to top), t=pings
% (time))
% x = vertical coordinate, in units of metres from the transducer (bottom
% to top)
% StartBin, Params, Model, and phi_increment are the same as inputs (see above)  
%
% Daniel Buscombe May-June 2012

for f=1:length(basepath)

files=ReadImDir(basepath{f},'mat');
addpath(basepath{f})

if isempty(usefiles{f})
    usefiles{f}=[1:size(files,1)];
end

    for i=1:length(usefiles{f})

        load(deblank(files(usefiles{f}(i),:)))

        x=data.rrcor(:,1); % distance from transducer in metres
        Size=zeros(size(data.Vrcor,1),length(data.Vrcor));
        Mass=zeros(size(data.Vrcor,1),length(data.Vrcor));
        
        PhiMin = -log (2000*Params{f}(3)) / log (2);
        PhiMax = -log (2000*Params{f}(2)) / log (2);
        Phi = PhiMax:-phi_increment:PhiMin;     %Generate Sediment sizes in fraction of Phi increments
        As = (2.^-Phi')/2000; %scale to radius in (m)

        if isnan(data.c)
            data.c=1.4975e+03;
        end
        
        if isnan(data.D)
            data.D=0.911;
        end
        
        % calculate form function, scattering cross section and sediment density
        [ff XX Density] = CalcFormFunction(Model,As,Params{f}(1),[data.Abs.Freq],data.c);
        
        h = waitbar(0,'Please wait...');
        for k=1:length(data.Vrcor) % k is columns ('time')
        
        % for debug
        %Freq=[data.Abs.Freq]; Kt=data.AbsKt; r=x; V=squeeze(data.Vrcor(:,k,:));
        
        % calculate size distribution and mass concentration distribution
        [Size(:,k) Mass(:,k)] = CalcSedimentSizeAndMass_efficient([data.Abs.Freq],data.AbsKt,x,squeeze(data.Vrcor(:,k,:)),...
            StartBin, ff, XX, Density, As); 
        
        waitbar(k/length(data.Vrcor),h)
        end
        close(h)

        Size = Size.* 2000;   %As Diameter in mm rather than in m
        
        save([basepath{f},filesep,deblank(files(usefiles{f}(i),:)),'_size_n_mass.mat'],'Size',...
            'Mass','x','StartBin','Params','Model','phi_increment')
        
    end

end


