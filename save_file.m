function [ ] = save_file( handles )
%SAVE_FILE Summary of this function goes here
%   Detailed explanation goes here
data = handles.data;
save(fullfile(handles.data.path, handles.data.file_name), 'data');
end

