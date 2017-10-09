function data = w2_plot_custom_erps_by_channel(file_name,condition_nr,stats,path_to_save,data)

%-------DESCRIPTION---------------
%ADD
%----------------------------------

%-------PATH MANAGEMENT-----------
%----------------------------------
%add path of preprocessing 
addpath(genpath(fullfile(data.ieeglab_path,'processing')));

%modifies paths to include parent directory
path_to_save = fullfile(data.parent_directory, path_to_save);

%create directory where the trimmed sets will be stored
if ~exist(path_to_save, 'dir')
  mkdir(path_to_save);
end

%assert if file exists???

%-------LOAD PARAMETERS------------
%----------------------------------
condition_nr = str2num(condition_nr);

%assert possible values are 1 or 2

%---------RUN---------------------
%---------------------------------
if condition_nr == 1
    plot_erps_custom_by_channel(file_name,stats,path_to_save);
else
    plot_erps_custom_2_mean_conditions_by_channel(file_name,stats,path_to_save)
end