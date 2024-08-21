clear all; close all; clc;
[path,~,~] = fileparts(mfilename('fullpath'));
workspacePath = strcat(path, filesep);
cd(workspacePath)
addpath(genpath(workspacePath));

fprintf('workspacePath: %s \n', workspacePath);

%files_list = {'D:\Persyst\Persyst_EEGs_24h_Laydat\Patient_1.lay'};
files_list = get_files_list_fr253();
files_list = {'D:\FrLayDatMergedByDay\FR_1096\1096_DEMO_DAY4_FFT_ONLY\1096_day21.lay'}

% files_list = {'D:\FrLayDatMergedByDay\FR_253\2004-03-19_25301102.lay'};
% files_list = {'D:\FrLayDatMergedByDay\FR_253\2004-03-20_25301102.lay'};
% files_list = {'D:\FrLayDatMergedByDay\FR_253\2004-03-21_25301102.lay'};
% files_list = {'D:\FrLayDatMergedByDay\FR_253\2004-03-22_25301102.lay'};
% files_list = {'D:\FrLayDatMergedByDay\FR_253\2004-03-23_25301102.lay'};

pat_name = 'FR1096_day21_EEA';

extra_files = 0;

[Summary,SleepStage]=SleepSEEG_dlp(files_list,extra_files);

stages_datetimes = datetime(SleepStage(:,2), 'ConvertFrom', 'datenum', 'Format','MM-dd-yy HH:mm:ss');
staging_results_cell = cat(2, cellstr(stages_datetimes), num2cell(SleepStage(:,3)));
staging_results_table = cell2table(staging_results_cell, "VariableNames",["DateTime" "SleepStage"]);

time_now_str = cellstr(datetime('now','TimeZone','local','Format','dd-MM-yyyy HH:mm:ss'));
time_now_str = strrep(time_now_str, '-', '_'); time_now_str = strrep(time_now_str, ':', '_'); time_now_str = strrep(time_now_str, ' ', '__'); time_now_str = time_now_str{1};


staging_csv_fn = strcat('Output\', pat_name, '_', 'MNI_ieegSleepStages_', time_now_str, '.csv');
delete(staging_csv_fn);
writetable(staging_results_table, staging_csv_fn,'Delimiter',',');

function files_list = get_files_list_fr253()
    data_path = "D:/FrLayDatOneHrMax/pat_FR_253/";
    files_list = {...
        "2004-03-26_25301102_0199.lay";...
        "2004-03-18_25301102_0000.lay";...
        "2004-03-18_25301102_0001.lay";...
        "2004-03-18_25301102_0002.lay";...
        "2004-03-18_25301102_0003.lay";...
        "2004-03-18_25301102_0004.lay";...
        "2004-03-18_25301102_0005.lay";...
        "2004-03-18_25301102_0006.lay";...
        "2004-03-18_25301102_0007.lay";...
        "2004-03-18_25301102_0008.lay";...
        "2004-03-18_25301102_0009.lay";...
        "2004-03-18_25301102_0010.lay";...
        "2004-03-18_25301102_0011.lay";...
        "2004-03-18_25301102_0012.lay";...
        "2004-03-19_25301102_0013.lay";...
        "2004-03-19_25301102_0014.lay";...
        "2004-03-19_25301102_0015.lay";...
        "2004-03-19_25301102_0016.lay";...
        "2004-03-19_25301102_0017.lay";...
        "2004-03-19_25301102_0018.lay";...
        "2004-03-19_25301102_0019.lay";...
        "2004-03-19_25301102_0020.lay";...
        "2004-03-19_25301102_0021.lay";...
        "2004-03-19_25301102_0022.lay";...
        "2004-03-19_25301102_0023.lay";...
        "2004-03-19_25301102_0024.lay";...
        "2004-03-19_25301102_0025.lay";...
        "2004-03-19_25301102_0026.lay";...
        "2004-03-19_25301102_0027.lay";...
        "2004-03-19_25301102_0028.lay";...
        "2004-03-19_25301102_0029.lay";...
        "2004-03-19_25301102_0030.lay";...
        "2004-03-19_25301102_0031.lay";...
        "2004-03-19_25301102_0032.lay";...
        "2004-03-19_25301102_0033.lay";...
        "2004-03-19_25301102_0034.lay";...
        "2004-03-19_25301102_0035.lay";...
        "2004-03-19_25301102_0036.lay";...
        "2004-03-19_25301102_0037.lay";...
        "2004-03-20_25301102_0038.lay";...
        "2004-03-20_25301102_0039.lay";...
        "2004-03-20_25301102_0040.lay";...
        "2004-03-20_25301102_0041.lay";...
        "2004-03-20_25301102_0042.lay";...
        "2004-03-20_25301102_0043.lay";...
        "2004-03-20_25301102_0044.lay";...
        "2004-03-20_25301102_0045.lay";...
        "2004-03-20_25301102_0046.lay";...
        "2004-03-20_25301102_0047.lay";...
        "2004-03-20_25301102_0048.lay";...
        "2004-03-20_25301102_0049.lay";...
        "2004-03-20_25301102_0050.lay";...
        "2004-03-20_25301102_0051.lay";...
        "2004-03-20_25301102_0052.lay";...
        "2004-03-20_25301102_0053.lay";...
        "2004-03-20_25301102_0054.lay";...
        "2004-03-20_25301102_0055.lay";...
        "2004-03-20_25301102_0056.lay";...
        "2004-03-20_25301102_0057.lay";...
        "2004-03-20_25301102_0058.lay";...
        "2004-03-20_25301102_0059.lay";...
        "2004-03-20_25301102_0060.lay";...
        "2004-03-20_25301102_0061.lay";...
        "2004-03-21_25301102_0062.lay";...
        "2004-03-21_25301102_0063.lay";...
        "2004-03-21_25301102_0064.lay";...
        "2004-03-21_25301102_0065.lay";...
        "2004-03-21_25301102_0066.lay";...
        "2004-03-21_25301102_0067.lay";...
        "2004-03-21_25301102_0068.lay";...
        "2004-03-21_25301102_0069.lay";...
        "2004-03-21_25301102_0070.lay";...
        "2004-03-21_25301102_0071.lay";...
        "2004-03-21_25301102_0072.lay";...
        "2004-03-21_25301102_0073.lay";...
        "2004-03-21_25301102_0074.lay";...
        "2004-03-21_25301102_0075.lay";...
        "2004-03-21_25301102_0076.lay";...
        "2004-03-21_25301102_0077.lay";...
        "2004-03-21_25301102_0078.lay";...
        "2004-03-21_25301102_0079.lay";...
        "2004-03-21_25301102_0080.lay";...
        "2004-03-21_25301102_0081.lay";...
        "2004-03-21_25301102_0082.lay";...
        "2004-03-21_25301102_0083.lay";...
        "2004-03-21_25301102_0084.lay";...
        "2004-03-21_25301102_0085.lay";...
        "2004-03-21_25301102_0086.lay";...
        "2004-03-22_25301102_0087.lay";...
        "2004-03-22_25301102_0088.lay";...
        "2004-03-22_25301102_0089.lay";...
        "2004-03-22_25301102_0090.lay";...
        "2004-03-22_25301102_0091.lay";...
        "2004-03-22_25301102_0092.lay";...
        "2004-03-22_25301102_0093.lay";...
        "2004-03-22_25301102_0094.lay";...
        "2004-03-22_25301102_0095.lay";...
        "2004-03-22_25301102_0096.lay";...
        "2004-03-22_25301102_0097.lay";...
        "2004-03-22_25301102_0098.lay";...
        "2004-03-22_25301102_0099.lay";...
        "2004-03-22_25301102_0100.lay";...
        "2004-03-22_25301102_0101.lay";...
        "2004-03-22_25301102_0102.lay";...
        "2004-03-22_25301102_0103.lay";...
        "2004-03-22_25301102_0104.lay";...
        "2004-03-22_25301102_0105.lay";...
        "2004-03-22_25301102_0106.lay";...
        "2004-03-22_25301102_0107.lay";...
        "2004-03-22_25301102_0108.lay";...
        "2004-03-22_25301102_0109.lay";...
        "2004-03-22_25301102_0110.lay";...
        "2004-03-23_25301102_0111.lay";...
        "2004-03-23_25301102_0112.lay";...
        "2004-03-23_25301102_0113.lay";...
        "2004-03-23_25301102_0114.lay";...
        "2004-03-23_25301102_0115.lay";...
        "2004-03-23_25301102_0116.lay";...
        "2004-03-23_25301102_0117.lay";...
        "2004-03-23_25301102_0118.lay";...
        "2004-03-23_25301102_0119.lay";...
        "2004-03-23_25301102_0120.lay";...
        "2004-03-23_25301102_0121.lay";...
        "2004-03-23_25301102_0122.lay";...
        "2004-03-23_25301102_0123.lay";...
        "2004-03-23_25301102_0124.lay";...
        "2004-03-23_25301102_0125.lay";...
        "2004-03-23_25301102_0126.lay";...
        "2004-03-23_25301102_0127.lay";...
        "2004-03-23_25301102_0128.lay";...
        "2004-03-23_25301102_0129.lay";...
        "2004-03-23_25301102_0130.lay";...
        "2004-03-23_25301102_0131.lay";...
        "2004-03-23_25301102_0132.lay";...
        "2004-03-23_25301102_0133.lay";...
        "2004-03-23_25301102_0134.lay";...
        "2004-03-24_25301102_0135.lay";...
        "2004-03-24_25301102_0136.lay";...
        "2004-03-24_25301102_0137.lay";...
        "2004-03-24_25301102_0138.lay";...
        "2004-03-24_25301102_0139.lay";...
        "2004-03-24_25301102_0140.lay";...
        "2004-03-24_25301102_0141.lay";...
        "2004-03-24_25301102_0142.lay";...
        "2004-03-24_25301102_0143.lay";...
        "2004-03-24_25301102_0144.lay";...
        "2004-03-24_25301102_0145.lay";...
        "2004-03-24_25301102_0146.lay";...
        "2004-03-24_25301102_0147.lay";...
        "2004-03-24_25301102_0148.lay";...
        "2004-03-24_25301102_0149.lay";...
        "2004-03-24_25301102_0150.lay";...
        "2004-03-24_25301102_0151.lay";...
        "2004-03-24_25301102_0152.lay";...
        "2004-03-24_25301102_0153.lay";...
        "2004-03-24_25301102_0154.lay";...
        "2004-03-24_25301102_0155.lay";...
        "2004-03-24_25301102_0156.lay";...
        "2004-03-24_25301102_0157.lay";...
        "2004-03-24_25301102_0158.lay";...
        "2004-03-25_25301102_0159.lay";...
        "2004-03-25_25301102_0160.lay";...
        "2004-03-25_25301102_0161.lay";...
        "2004-03-25_25301102_0162.lay";...
        "2004-03-25_25301102_0163.lay";...
        "2004-03-25_25301102_0164.lay";...
        "2004-03-25_25301102_0165.lay";...
        "2004-03-25_25301102_0166.lay";...
        "2004-03-25_25301102_0167.lay";...
        "2004-03-25_25301102_0168.lay";...
        "2004-03-25_25301102_0169.lay";...
        "2004-03-25_25301102_0170.lay";...
        "2004-03-25_25301102_0171.lay";...
        "2004-03-25_25301102_0172.lay";...
        "2004-03-25_25301102_0173.lay";...
        "2004-03-25_25301102_0174.lay";...
        "2004-03-25_25301102_0175.lay";...
        "2004-03-25_25301102_0176.lay";...
        "2004-03-25_25301102_0177.lay";...
        "2004-03-25_25301102_0178.lay";...
        "2004-03-25_25301102_0179.lay";...
        "2004-03-25_25301102_0180.lay";...
        "2004-03-25_25301102_0181.lay";...
        "2004-03-25_25301102_0182.lay";...
        "2004-03-26_25301102_0183.lay";...
        "2004-03-26_25301102_0184.lay";...
        "2004-03-26_25301102_0185.lay";...
        "2004-03-26_25301102_0186.lay";...
        "2004-03-26_25301102_0187.lay";...
        "2004-03-26_25301102_0188.lay";...
        "2004-03-26_25301102_0189.lay";...
        "2004-03-26_25301102_0190.lay";...
        "2004-03-26_25301102_0191.lay";...
        "2004-03-26_25301102_0192.lay";...
        "2004-03-26_25301102_0193.lay";...
        "2004-03-26_25301102_0194.lay";...
        "2004-03-26_25301102_0195.lay";...
        "2004-03-26_25301102_0196.lay";...
        "2004-03-26_25301102_0197.lay";...
        "2004-03-26_25301102_0198.lay";...
    };
    %files_list = cellfun(@(c)strcat(data_path,c), files_list, 'uni', false);
    files_list = cellfun(@(c)strcat(data_path,c), files_list, 'UniformOutput', false);
    %files_list = cellfun(@(c)strcat(data_path,c), files_list, 'UniformOutput', false);

    filesList = cellfun(@(x) x{1}, files_list, 'UniformOutput', false);
    filesList = sort(filesList);
    files_list = filesList;
end