function [EEG,data] = w_erps_by_channel_EEG(path_to_files,condition_1,condition_2,cycles,freq_range,alpha,fdr,scale,basenorm,erps_max,path_to_save,EEG,data)

data.parent_directory = '';
[data] = w2_erps_by_channel_EEG(path_to_files,condition_1,condition_2,cycles,freq_range,alpha,fdr,scale,basenorm,erps_max,path_to_save,data);