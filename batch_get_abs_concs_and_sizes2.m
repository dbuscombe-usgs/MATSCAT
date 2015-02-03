
% batch_get_abs_concs_and_sizes2
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
basepath{1}='/home/hercules/Desktop/nov8_withfobs/data_out/processed';
% basepath{2}='/media/CRUZER/repeat_abs_fobs_5nov/data_out/processed';


usefiles{1}=[]; %[]; % index of files 
usefiles{2}=[];
% usefiles{3}=[]; % leave blank if you want to process all the files in the stated directory
% usefiles{4}=[];
% usefiles{5}=[];
% usefiles{6}=[];

% RADIUS NOT DIAMETER
% StDev, min size in metres, max size in metres
Params{1}=[0.15    0.0002/2   0.0004/2]; %[0.15    0.0001/2    0.00015/2]; 
Params{2}=[0.15    0.0002/2   0.0004/2]; %[0.15    0.0001/2    0.00015/2]; 
% Params{2}=[0.15    0.00011/2   0.0006/2];
% Params{3}=[0.15    0.00014/2    0.00021/2];
% Params{4}=[0.15    0.00014/2    0.00021/2];
% 
% Params{5}=[0.15    0.000063/2    0.00015/2];
% Params{6}=[0.15    0.000063/2    0.00015/2];


StartBin=10; % start calcs at bin 10
phi_increment=.05; %resolution on sediment distribution in phi units
Model='GLASS'; % 'GLASS' for glass spheres, 'SAND' for sand

% pass these arguments to the function
get_abs_concs_and_sizes(basepath,usefiles,Params,StartBin,phi_increment,Model)

