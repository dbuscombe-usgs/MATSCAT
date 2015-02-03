
% batch_range_correct_abs
close all; clear all;clc

% example wrapper function to call range_correct_abs.m in batch mode
% Daniel Buscombe May-June 2012
% version 1.0

if isunix
    if ~exist('CalcSedimentSizeAndMass_efficient.m','file')
        [stat,path]=system('locate CalcSedimentSizeAndMass_efficient.m');
        addpath(fileparts(path))
    end
else
    if ~exist('CalcSedimentSizeAndMass_efficient.m','file')    
        % WINDOWS USERS EDIT THIS TO WHERE THE ABS TOOLBOX IS
        % IF CalcSedimentSizeAndMass_efficient.m is not in pwd
        addpath(path_to_the_abs_toolbox)    
    end
end

% directories to process
basepath{1}='/home/hercules/Desktop/nov8_withfobs';
% basepath{2}='/media/CRUZER/repeat_abs_fobs_5nov';
%basepath{2}='~/Desktop/lab_data_aug2012_raw/abs/aug20';
%basepath{3}='~/Desktop/lab_data_aug2012_raw/abs/aug23';
%basepath{4}='~/Desktop/lab_data_aug2012_raw/abs/aug29';


temp{1}=12;%Temperature in degrees C
depth{1}=0.6; % depth in m
salinity{1}=0; % salinity

% attenuation coefficients as determined by manufacturer
AbsKt{1}(1) = 0.01288;
AbsKt{1}(2) = 0.01747;
AbsKt{1}(3) = 0.00176;

temp{2}=temp{1};
depth{2}=depth{1};
salinity{2}=salinity{1};
AbsKt{2}=AbsKt{1};

temp{3}=temp{1};
depth{3}=depth{1};
salinity{3}=salinity{1};
AbsKt{3}=AbsKt{1};

temp{4}=temp{1};
depth{4}=depth{1};
salinity{4}=salinity{1};
AbsKt{4}=AbsKt{1};


% pass arguments to function
range_correct_abs(basepath,temp,depth,salinity,AbsKt)


