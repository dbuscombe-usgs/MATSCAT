function [filenames] = ReadImDir(directory,ext)

%READIMDIR 
% Reads images in a directory, and returns them as a structure of filenames
% as strings
% [filenames] = ReadImDir(directory,extension)

% written by Daniel Buscombe, 07/2006
% version 1.0

direc=dir([directory,filesep,'*.',ext]); %list directory and separate .*ext files
filenames={};   %create a structure of these files
[filenames{1:length(direc),1}] = deal(direc.name); %Deal inputs to outputs!
filenames=sortrows(char(filenames{:})); %Create character array, and sort rows in ascending order

return
 
%sometimes needs to be inputted into 'sortn.m' to make sure numbers in correct order when there are
%also letters to confuse things
