clear all; close all; clc;
[path,~,~] = fileparts(mfilename('fullpath'));
workspacePath = strcat(path, filesep);
cd(workspacePath)
addpath(genpath(workspacePath));

fprintf('workspacePath: %s \n', workspacePath);

FileList = {"C:\Users\DLP\Documents\Elpi_Exports\JOHNSON~ Isla_c9e2488f-5af8-4947-a4cc-53f7b99a90d9_Exported_2024-7-29_19-8-12.edf"};
FileList = {'D:\Persyst\Persyst_EEGs_24h_Laydat\Patient_1.lay'};
ExtraFiles = 0;

[Summary,SleepStage]=SleepSEEG_dlp(FileList,ExtraFiles);