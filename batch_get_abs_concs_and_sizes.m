
% batch_get_abs_concs_and_sizes
close all; clear all;clc

% example wrapper function to call get_abs_concs_and_sizes.m in batch mode
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

% specify a number of directories to process
basepath{1}='~/Desktop/lab_data_aug2012_raw/abs/aug1/data_out/processed';
% basepath{2}='~/Desktop/lab_data_aug2012_raw/abs/aug20/data_out/processed';
% basepath{3}='~/Desktop/lab_data_aug2012_raw/abs/aug23/data_out/processed';
% basepath{4}='~/Desktop/lab_data_aug2012_raw/abs/aug29/data_out/processed';


usefiles{1}=[]; %[1:7,11]; % index of files 
% usefiles{2}=[];
% usefiles{3}=[]; % leave blank if you want to process all the files in the stated directory
% usefiles{4}=[];
% usefiles{5}=[];
% usefiles{6}=[];

% RADIUS NOT DIAMETER
% StDev, min size in metres, max size in metres
Params{1}=[0.15    0.00001/2   0.0002/2]; %[0.15    0.0001/2    0.00015/2]; 
% Params{2}=[0.15    0.00001/2   0.0002/2];
% Params{3}=[0.15    0.00001/2    0.0002/2];
% Params{4}=[0.15    0.00001/2    0.0002/2];
% 
% Params{5}=[0.15    0.000063/2    0.00015/2];
% Params{6}=[0.15    0.000063/2    0.00015/2];


StartBin=10; % start calcs at bin 10
phi_increment=.05; %resolution on sediment distribution in phi units
Model='SAND'; % 'GLASS' and 'SAND' for glass spheres

% pass these arguments to the function
get_abs_concs_and_sizes(basepath,usefiles,Params,StartBin,phi_increment,Model)

