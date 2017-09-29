function [ handles ] = load_file( file_name, path, handles )
%SAVE_FILE Summary of this function goes here
%   Detailed explanation goes here
load(fullfile(path, file_name));
[pathstr,~,~] = fileparts(which('ieeglab'));
data.ieeglab_path = pathstr;
handles.data = data;
end

