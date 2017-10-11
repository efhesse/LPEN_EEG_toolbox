function [EEG,data] = w_custom_newtimef_2_precalculated_conditions_by_channel_EEG(path_to_files,condition1_data,condition2_data,path_to_save,EEG,data)

data.parent_directory = '';
[data] = w2_custom_newtimef_2_precalculated_conditions_by_channel_EEG(path_to_files,condition1_data,condition2_data,path_to_save,data);