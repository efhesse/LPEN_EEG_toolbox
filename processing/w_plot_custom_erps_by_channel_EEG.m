function [EEG,data] = w_plot_custom_erps_by_channel_EEG(path_to_file,file_name,condition_nr,stats,path_to_save,EEG,data)

data.parent_directory = '';
[data] = w2_plot_custom_erps_by_channel_EEG(path_to_file,file_name,condition_nr,stats,path_to_save,data);